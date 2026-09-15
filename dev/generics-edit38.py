# N8: US English. `summarise`, `summariser` and `summarises` are
# British. The historical dev/generics-edit*.R scripts are left as
# they are, because each one records the exact text it replaced.
# "Nature Human Behaviour" in NEWS.md is a journal's name and is not
# this lane's.
import io

ROOT = "C:/Users/adf44/source/r/frmtmb-wt-generics/"
edits = {
    "R/scales.R": [("summarises", "summarizes")],
    "dev/generics-findings.md": [("summariser", "summarizer"),
                                 ("summarises", "summarizes"),
                                 ("summarise ", "summarize ")],
    "dev/generics-assemble.R": [("summariser", "summarizer")],
    "dev/generics-summary.R": [("summariser", "summarizer")],
}
for rel, pairs in edits.items():
    p = ROOT + rel
    s = io.open(p, encoding="utf-8", newline="").read()
    for old, new in pairs:
        n = s.count(old)
        if n == 0:
            raise SystemExit("no %r in %s" % (old, rel))
        s = s.replace(old, new)
        print("%s: %d x %s -> %s" % (rel, n, old, new))
    io.open(p, "w", encoding="utf-8", newline="\n").write(s)
