# Item 764 widened: every local_mocked_bindings() call in the test files
# of core and the extensions that names no .package (comments skipped
# inside a call; a match inside a comment line is a false positive).
#   uv run --no-project python dev/hygiene-mockscan.py
import re, glob, sys
files = glob.glob("tests/testthat/test-*.R") + glob.glob("extensions/*/tests/testthat/test-*.R")
tot = 0; per = {}
for f in sorted(files):
    s = open(f, encoding="utf-8").read()
    for m in re.finditer(r"local_mocked_bindings\(", s):
        # find the matching close paren
        i = m.end(); depth = 1; q = None
        while i < len(s) and depth:
            c = s[i]
            if q:
                if c == "\\": i += 1
                elif c == q: q = None
            elif c in "\"'": q = c
            elif c == "#":
                j = s.find("\n", i); i = j if j > 0 else len(s); continue
            elif c == "(": depth += 1
            elif c == ")": depth -= 1
            i += 1
        call = s[m.start():i]
        line = s.count("\n", 0, m.start()) + 1
        if ".package" not in call:
            tot += 1; per.setdefault(f, []).append(line)
for f, ls in per.items(): print(f, ls)
print("calls without .package:", tot, "in", len(per), "files")
