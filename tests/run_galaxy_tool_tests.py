#!/usr/bin/env python
"""Run a Galaxy tool's <tests> without a Galaxy server.

`planemo test` is the reference check and CI runs it, but it needs to start a
Galaxy instance on a local port, which is not possible in every development
sandbox. This runner covers the part of the test that concerns the wrapper:
it renders the tool's <command> with Cheetah from the test's parameter values
plus the XML defaults for everything the test leaves out, executes the
rendered command in a scratch work directory, and checks the declared
<assert_contents> and expect_num_outputs against the files the tool wrote.

What it does NOT check: datatype sniffing, metadata, output format
declarations, dependency resolution, or the diff-style output comparisons.
Those need the real Galaxy job runner - use planemo for them.

Usage:
    python tests/run_galaxy_tool_tests.py galaxy/ewas_dmr_ml.xml [...]
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile

from Cheetah.Template import Template
from galaxy.tool_util.parser import get_tool_source


class Param(str):
    """A parameter that renders as its value and carries child attributes."""

    def __new__(cls, value, children=None, truth=None):
        o = super().__new__(cls, value)
        o._children = children or {}
        o._truth = truth
        return o

    def __getattr__(self, name):
        try:
            return self._children[name]
        except KeyError:
            raise AttributeError(name)

    def __bool__(self):
        return bool(self) if self._truth is None else self._truth

    # Cheetah calls this for $x usage
    def __str__(self):
        return str.__str__(self)


def _bool_param(elem):
    checked = (elem.get("checked", "false") or "").lower() in ("true", "yes", "1")
    truevalue = elem.get("truevalue", "true")
    falsevalue = elem.get("falsevalue", "false")
    return Param(truevalue if checked else falsevalue, truth=checked)


def _select_default(elem):
    opts = elem.findall("option")
    for o in opts:
        if (o.get("selected", "false") or "").lower() == "true":
            return Param(o.get("value", ""))
    return Param(opts[0].get("value", "")) if opts else Param("")


def collect_defaults(elem):
    """Defaults for one <inputs>/<section> subtree, as a nested dict."""
    out = {}
    for child in elem:
        name = child.get("name")
        if child.tag == "section":
            out[name] = Param("", children=collect_defaults(child))
        elif child.tag == "conditional":
            sel = child.find("param")
            chosen = _select_default(sel) if sel is not None else Param("")
            kids = {sel.get("name"): chosen} if sel is not None else {}
            for when in child.findall("when"):
                if when.get("value") == str(chosen):
                    kids.update(collect_defaults(when))
            out[name] = Param(str(chosen), children=kids)
        elif child.tag == "repeat":
            out[name] = []
        elif child.tag == "param":
            ptype = child.get("type", "text")
            if ptype == "data":
                continue  # only a test can supply one
            elif ptype == "boolean":
                out[name] = _bool_param(child)
            elif ptype == "select":
                out[name] = _select_default(child)
            else:
                out[name] = Param(child.get("value", ""))
    return out


def apply_test_values(ns, inputs, data_dir, boolean_names):
    """Overlay a test's flat `a|b` parameter values onto the defaults."""
    for item in inputs:
        parts = item["name"].split("|")
        value = item["value"]
        if value is None:
            value = ""
        # a data param's value is a file name under test-data/
        candidate = os.path.join(data_dir, str(value))
        if os.path.exists(candidate):
            value = os.path.abspath(candidate)
        target = ns
        for p in parts[:-1]:
            target = target[p]._children if isinstance(target, dict) else target._children
        leaf = parts[-1]
        if "|".join(parts) in boolean_names or leaf in boolean_names:
            truth = str(value).lower() in ("true", "yes", "1")
            target[leaf] = Param(str(value), truth=truth)
        else:
            target[leaf] = Param(str(value))
    return ns


def boolean_param_names(elem, prefix=""):
    names = set()
    for child in elem:
        name = child.get("name") or ""
        full = f"{prefix}{name}"
        if child.tag in ("section", "conditional", "when", "repeat"):
            names |= boolean_param_names(child, f"{full}|" if name else prefix)
        elif child.tag == "param" and child.get("type") == "boolean":
            names.add(full)
            names.add(name)
    return names


def check_assertions(path, assert_list):
    """Only the assertion tags these wrappers use."""
    problems = []
    try:
        text = open(path, errors="replace").read()
    except OSError as e:
        return [f"unreadable: {e}"]
    for a in assert_list:
        tag, attrs = a["tag"], a["attributes"]
        if tag == "has_text":
            if attrs["text"] not in text:
                problems.append(f"has_text {attrs['text']!r} absent")
        elif tag == "has_text_matching":
            if not re.search(attrs["expression"], text):
                problems.append(f"has_text_matching {attrs['expression']!r} no match")
        elif tag == "has_n_lines":
            n = len(text.splitlines())
            if int(attrs.get("n", -1)) not in (-1, n):
                problems.append(f"has_n_lines {attrs['n']} != {n}")
        else:
            problems.append(f"unsupported assertion <{tag}> - use planemo for this test")
    return problems


def run_tool_tests(tool_xml, keep=False):
    ts = get_tool_source(tool_xml)
    root = ts.root
    tool_dir = os.path.abspath(os.path.dirname(tool_xml))
    data_dir = os.path.join(tool_dir, "test-data")
    inputs_elem = root.find("inputs")
    bools = boolean_param_names(inputs_elem)
    outputs = ts.parse_outputs(None)[0]
    command_tpl = ts.parse_command()
    tests = ts.parse_tests_to_dict()["tests"]
    tool_id = ts.parse_id()

    if not tests:
        print(f"{tool_id}: no <tests> defined - nothing to run")
        return 0

    failures = 0
    for i, t in enumerate(tests, 1):
        ns = collect_defaults(inputs_elem)
        ns = apply_test_values(ns, t["inputs"], data_dir, bools)
        ns["__tool_directory__"] = Param(tool_dir)
        work = tempfile.mkdtemp(prefix=f"{tool_id}_t{i}_")
        # Galaxy collapses newlines in the rendered command to spaces before
        # handing it to the shell, which is why wrappers may put one argument
        # per line. Do the same, or every continuation line runs as a command.
        rendered = str(Template(source=command_tpl, searchList=[ns]))
        rendered = re.sub(r"[\r\n]+", " ", rendered).strip()
        with open(os.path.join(work, "command.sh"), "w") as fh:
            fh.write(rendered)
        proc = subprocess.run(["bash", "command.sh"], cwd=work,
                              capture_output=True, text=True)
        problems = []
        if proc.returncode != int(t.get("expect_exit_code") or 0):
            problems.append(f"exit code {proc.returncode}")
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-8:]
            problems += [f"  | {line}" for line in tail]

        produced = 0
        for name, out in outputs.items():
            fwd = getattr(out, "from_work_dir", None)
            path = os.path.join(work, fwd) if fwd else os.path.join(work, name)
            if os.path.exists(path) and os.path.getsize(path) > 0:
                produced += 1
            else:
                problems.append(f"output {name}: missing or empty ({fwd})")
        expect_n = t.get("expect_num_outputs")
        if expect_n is not None and produced != int(expect_n):
            problems.append(f"expect_num_outputs {expect_n}, produced {produced}")

        for o in t["outputs"]:
            out = outputs[o["name"]]
            fwd = getattr(out, "from_work_dir", None)
            path = os.path.join(work, fwd) if fwd else os.path.join(work, o["name"])
            problems += [f"output {o['name']}: {p}"
                         for p in check_assertions(path, o["attributes"]["assert_list"])]

        if problems:
            failures += 1
            print(f"FAIL {tool_id} test {i}  (work dir kept: {work})")
            for p in problems:
                print("   ", p)
        else:
            print(f"ok   {tool_id} test {i}")
            if not keep:
                shutil.rmtree(work, ignore_errors=True)
    return failures


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    keep = "--keep" in sys.argv[1:]
    sys.exit(min(1, sum(run_tool_tests(a, keep) for a in args)))
