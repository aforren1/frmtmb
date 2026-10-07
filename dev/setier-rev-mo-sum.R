# Reviewer of lane setier: compare the 200-seed mo() study, lane vs base.
#   Rscript dev/setier-rev-mo-sum.R
rd <- function(arm) {
  do.call(rbind, lapply(1:4, function(c) {
    read.delim(sprintf("dev/setier-rev-log/mo/mo-%s-%d.tsv", arm, c),
               stringsAsFactors = FALSE, quote = "")
  }))
}
a <- rd("base")
b <- rd("lane")
stopifnot(identical(a$seed, b$seed), identical(a$form, b$form))
cat("rows", nrow(a), "seeds", length(unique(a$seed)), "\n")
for (x in list(list("base", a), list("lane", b))) {
  d <- x[[2]]
  lostfits <- d$n_lost > 0 | d$se_finite < d$n_par
  cat(x[[1]], ": lost-SE fits", sum(lostfits),
      " SE warning", sum(d$se_warn & lostfits),
      " code!=0", sum(d$code != 0 & lostfits),
      " silent", sum(lostfits & !d$se_warn & d$code == 0 &
                       d$n_warn_fit == 0),
      " false alarms (warn, nothing lost)", sum(d$se_warn & !lostfits),
      " ce fail", sum(!d$ce_ok),
      " max rel b (finite)", signif(max(d$max_rel_b[is.finite(d$max_rel_b)]),
                                    3), "\n")
}
cols <- setdiff(names(a), c("se_b", "max_rel_b", "se_b_ref"))
diff <- which(Reduce(`|`, lapply(cols, function(c) {
  !(a[[c]] == b[[c]] | (is.na(a[[c]]) & is.na(b[[c]])))
})))
cat("rows differing in", paste(cols, collapse = ","), ":", length(diff), "\n")
for (i in diff) {
  cat(sprintf("seed %d %s | base: lost=%s warn=%s | lane: lost=%s warn=%s\n",
              a$seed[i], a$form[i], a$lost[i], substr(a$warn_fit[i], 1, 90),
              b$lost[i], substr(b$warn_fit[i], 1, 90)))
}
nb <- sum(a$se_b != b$se_b)
cat("rows whose coefficient SEs differ:", nb, "\n")
