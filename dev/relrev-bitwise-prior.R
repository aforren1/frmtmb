# Reviewer: the 30 default_prior() outputs that differ between r5 and r6
# in dev/relrev-bitwise.R: are they exactly the added per-threshold
# Intercept rows (lane fixes)? And did the one-step residuals run?
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev-log"
a <- readRDS(file.path(root, "relrev-bitwise-r5.rds"))
b <- readRDS(file.path(root, "relrev-bitwise-r6.rds"))
for (m in names(a)) {
  pa <- as.data.frame(a[[m]]$prior); pb <- as.data.frame(b[[m]]$prior)
  if (identical(a[[m]]$prior, b[[m]]$prior)) next
  key <- function(p) do.call(paste, c(p[, intersect(c("prior", "class", "coef", "group",
                                                     "resp", "dpar", "nlpar", "lb", "ub",
                                                     "source"), names(p))], sep = "|"))
  ka <- key(pa); kb <- key(pb)
  added <- pb[!kb %in% ka, ]
  lost <- pa[!ka %in% kb, ]
  only_thr <- all(added$class == "Intercept" & nzchar(added$coef))
  cat(sprintf("%-32s rows %d -> %d; added %d (all class Intercept with coef: %s); lost %d; cols same: %s\n",
              m, nrow(pa), nrow(pb), nrow(added), only_thr, nrow(lost),
              identical(names(pa), names(pb))))
}
osa <- vapply(names(a), function(m) {
  if (is.null(a[[m]]$osa)) return(NA_character_)
  if (is.character(a[[m]]$osa)) return(paste("ERR", a[[m]]$osa))
  paste(identical(a[[m]]$osa, b[[m]]$osa), "finite", all(is.finite(as.matrix(a[[m]]$osa[, 1]))))
}, "")
print(osa[!is.na(osa)])
