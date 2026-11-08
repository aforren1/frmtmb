"""Reviewer check of item 3: parse every workflow, list each step that
runs skip-gate.sh with its `if`, the id it reads, and whether a step of
that id precedes it in the same job with no `if` of its own and no
continue-on-error. Also job and step timeouts.

    uv run --no-project --with pyyaml python dev/review-ci-yaml.py <dir>
"""
import pathlib
import re
import sys

import yaml

root = pathlib.Path(sys.argv[1])
n_gates = 0
for p in sorted(root.glob(".github/workflows/*.yaml")):
    raw = p.read_bytes()
    crlf = raw.count(b"\r\n")
    lf = raw.count(b"\n")
    doc = yaml.safe_load(raw)
    for jname, job in doc["jobs"].items():
        steps = job.get("steps", [])
        jt = job.get("timeout-minutes")
        for i, s in enumerate(steps):
            run = s.get("run", "") or ""
            if "skip-gate.sh" not in run:
                continue
            n_gates += 1
            cond = s.get("if")
            m = re.search(r"steps\.(\w+)\.outcome", cond or "")
            sid = m.group(1) if m else None
            prior = [t for t in steps[:i] if t.get("id") == sid]
            tgt = prior[0] if prior else {}
            between = [t.get("name") or t.get("uses") for t in steps[steps.index(tgt) + 1:i]] if prior else []
            print(f"{p.name}:{jname}: if={cond!r}")
            print(f"    target id={sid} found={bool(prior)} "
                  f"uses={tgt.get('uses') or tgt.get('name')} "
                  f"target_if={tgt.get('if')!r} "
                  f"coe={tgt.get('continue-on-error')!r} "
                  f"step_timeout={tgt.get('timeout-minutes')} job_timeout={jt}")
            print(f"    steps between: {between}")
    print(f"  [{p.name}] parsed; CRLF {crlf} of {lf} line ends")
print("gates:", n_gates)
