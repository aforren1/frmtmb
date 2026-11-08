"""List the test files of one package that can reach stan_cache_write():
helper functions that call it, transitively, then the test files that
call any of them.

    uv run --no-project python -I dev/ci-cachecallers.py <tests/testthat dir>
"""
import pathlib
import re
import sys

d = pathlib.Path(sys.argv[1])
defn = re.compile(r"^([A-Za-z_.][A-Za-z0-9_.]*) <- function", re.M)
bodies = {}
for h in sorted(d.glob("*.R")):
    text = h.read_text(encoding="utf-8", errors="replace")
    starts = [(m.start(), m.group(1)) for m in defn.finditer(text)]
    for i, (pos, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        bodies[(h.name, name)] = text[pos:end]
names = {"stan_cache_write"}
while True:
    grown = {n for (f, n), b in bodies.items()
             if f.startswith("helper") and n not in names
             and any(re.search(r"\b" + re.escape(x) + r"\(", b) for x in names)}
    if not grown:
        break
    names |= grown
pat = re.compile(r"\b(" + "|".join(re.escape(n) for n in names) + r")\(")
files = sorted(t.name for t in d.glob("test-*.R")
               if pat.search(t.read_text(encoding="utf-8", errors="replace")))
print("functions:", " ".join(sorted(names)))
print(len(files), "files")
for f in files:
    print(f)
