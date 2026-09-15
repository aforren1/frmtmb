# Run after dev/generics-edit43.py. The findings page was assembled in
# round 2, so two of its generated blocks sit inline with round-2
# numbers. Turn those fences back into markers so dev/generics-assemble.R
# splices the round-3 blocks instead of leaving stale counts on the page.
import re

P = "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-findings.md"
s = open(P, encoding="utf-8", newline="").read()


def refence(title_line, marker):
    global s
    pat = re.compile(r"```\n" + re.escape(title_line) + r".*?\n```\n", re.S)
    hits = pat.findall(s)
    if len(hits) != 1:
        raise SystemExit("block %r found %d times" % (title_line, len(hits)))
    s = pat.sub("<!--BLOCK:" + marker + "-->\n", s, count=1)
    print("refenced", marker)


refence("== suite: frmtmb.sample, 17 files, one process per file, "
        "against the final core ==", "suitesample")
refence("== R CMD check --as-cran, built WITH vignettes, no --no-manual ==",
        "check")
open(P, "w", encoding="utf-8", newline="\n").write(s)
