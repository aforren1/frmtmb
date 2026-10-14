# roxygenise frmtmb.eam in the worktree only (man/ and NAMESPACE), and
# render ?wiener to text to check the new sections.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
p <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam"
roxygen2::roxygenise(p)
tools::Rd2txt(file.path(p, "man", "wiener.Rd"),
  out = "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out/wiener-rd.txt")
