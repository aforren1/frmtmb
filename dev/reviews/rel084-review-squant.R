# rel084 release review, squant fix round 2: the 12 group-level compat rows
# of quantile_elf are "conditional" with the note that points at
# ?quantile_elf "Group-level terms", in core alone and with every extension
# loaded; the totals (84 rows: works 26, conditional 44, refused 13,
# untested 1); the help section exists and says smoothing = 1.
#   REVLIB=r23 Rscript dev/reviews/rel084-review-squant.R [all]
source("dev/reviews/rel084-review-common.R")
if (identical(commandArgs(TRUE)[1], "all")) {
  for (p in c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
              "frmtmb.ode", "frmtmb.sample", "frmtmb.spatial", "frmtmb.survival")) {
    suppressMessages(library(p, character.only = TRUE))
  }
  cat("loaded:", paste(sort(grep("^frmtmb", loadedNamespaces(), value = TRUE)),
                       collapse = " "), "\n")
}
q <- frm_compat("quantile_elf")
o <- ifelse(q$feature_a == "quantile_elf", q$feature_b, q$feature_a)
cat("rows:", nrow(q), "\n")
print(table(q$status))
named <- c("us", "diag", "homdiag", "cs", "toep", "homtoep", "homcs", "us_t",
           "diag_t", "s()", "t2()", "gp_pred()")
for (nm in named) {
  i <- which(o == nm)
  cat(sprintf("%-10s %s %s | %s\n", nm, length(i),
              paste(q$status[i], collapse = ","),
              substr(q$note[i[1]], 1, 70)))
}
cat("all 12 conditional:", all(vapply(named, function(nm)
  all(q$status[o == nm] == "conditional"), NA)), "\n")
# the analogous asym_laplace rows, for comparison
a <- frm_compat("asym_laplace")
oa <- ifelse(a$feature_a == "asym_laplace", a$feature_b, a$feature_a)
cat("asym_laplace:", paste(sprintf("%s=%s", named, vapply(named, function(nm)
  paste(a$status[oa == nm], collapse = ","), "")), collapse = " "), "\n")
rd <- tools::Rd_db("frmtmb")
qrd <- rd[[grep("quantile_elf", names(rd))[1]]]
txt <- paste(capture.output(tools::Rd2txt(qrd, options = list(underline_titles = FALSE))),
             collapse = "\n")
cat("Rd file:", grep("quantile_elf", names(rd), value = TRUE)[1], "\n")
cat("has 'Group-level terms':", grepl("Group-level terms", txt), "\n")
seg <- regmatches(txt, regexpr("Group-level terms[\\s\\S]{0,1800}", txt, perl = TRUE))
cat(seg, "\n")
cat("DONE\n")
