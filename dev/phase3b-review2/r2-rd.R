# Reviewer 2, item 10: render the lane's Rd pages and grep the output.
for (f in c("extensions/frmtmb.eam/man/wiener.Rd", "extensions/frmtmb.eam/man/rdm.Rd",
            "extensions/frmtmb.eam/man/lba.Rd", "extensions/frmtmb.learn/man/bandit2arm_delta.Rd")) {
  out <- file.path("dev/phase3b-review2", paste0("rd-", basename(f), ".txt"))
  tools::Rd2txt(f, out = out, options = list(underline_titles = FALSE))
  cat(f, "->", out, length(readLines(out)), "lines\n")
}
