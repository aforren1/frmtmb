# The previous patch to dev/generics-assemble.R went through a bash
# heredoc, which ate a backslash and put literal newlines inside two R
# string literals. It parsed and would have worked, but it is not the
# code anyone would read. Written from a file this time.
P = "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-assemble.R"
s = open(P, encoding="utf-8").read()
bad = 'writeBin(charToRaw(paste0(paste(out, collapse = "\n"), "\n")), md)'
good = 'writeBin(charToRaw(paste0(paste(out, collapse = "\\n"), "\\n")), md)'
assert s.count(bad) == 1, "pattern not found"
s = s.replace(bad, good)
open(P, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
