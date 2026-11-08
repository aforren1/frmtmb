"""Make every skip-gate step run on a failed check too (lane ci, item 3).

Gives the step that writes testthat.Rout an id and conditions the gate
on that step having run, whatever its result. Run once from the
worktree root; it refuses a file it has already edited.
"""
import pathlib
import sys

GATE = ("      - name: Fail on a skip from a broken or missing Suggests"
        " package\n")
COND = (
    "      # on a failed check too: the gate ran only after a pass,"
    " so 0.79.0's\n"
    "      # other failures hid a skip on an extension until 0.79.1\n"
)


def gate_if(step_id):
    return (
        "        if: ${{ !cancelled() && steps." + step_id
        + ".outcome != 'skipped' }}\n"
    )


def edit(path, anchor, step_id, ngates):
    text = path.read_text(encoding="utf-8")
    if "steps." + step_id + ".outcome != 'skipped'" in text:
        sys.exit(f"{path}: already edited")
    lines = text.splitlines(keepends=True)
    out = []
    anchors = gates = 0
    for ln in lines:
        if ln == GATE:
            out.append(COND)
            out.append(ln)
            out.append(gate_if(step_id))
            gates += 1
            continue
        out.append(ln)
        if ln.strip() == anchor:
            indent = ln[: len(ln) - len(ln.lstrip())]
            # the id belongs to the step's mapping, two columns in from "-"
            out.append(indent + "  id: " + step_id + "\n")
            anchors += 1
    if gates != ngates or anchors != ngates:
        sys.exit(f"{path}: {gates} gates, {anchors} anchors, want {ngates}")
    path.write_text("".join(out), encoding="utf-8", newline="")
    print(f"{path}: {gates} gate(s)")


wf = pathlib.Path(".github/workflows")
for name in ["R-CMD-check.yaml"] + [
    f"check-frmtmb-{p}.yaml"
    for p in ["coupling", "eam", "latent", "learn", "ode", "sample",
              "spatial", "survival"]
]:
    edit(wf / name, "- uses: r-lib/actions/check-r-package@v2", "check", 1)
edit(wf / "test-coverage.yaml", "- name: Measure coverage", "measure", 2)
