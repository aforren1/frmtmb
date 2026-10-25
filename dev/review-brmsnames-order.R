## Prior-table row order, brms vs frmtmb (lane), from the saved outputs of
## dev/review-brmsnames-brms.R and dev/review-brmsnames-frm.R.
for (i in c(1:9, 13:15)) {
  b <- readRDS(sprintf("dev/review-brmsnames-out/brms-%02d.rds", i))
  f <- readRDS(sprintf("dev/review-brmsnames-out/frm-lane-%02d.rds", i))
  kb <- with(b$prior, paste(class, coef, group, resp, dpar, nlpar, sep = "|"))
  kf <- with(f$prior, paste(class, coef, group, resp, dpar, nlpar, sep = "|"))
  kf <- kf[kf %in% kb]
  kb <- kb[kb %in% kf]
  cat("== case", i, b$id, "same order on shared rows:", identical(kb, kf),
      "\n")
  if (!identical(kb, kf)) {
    print(cbind(brms = kb, frm = kf))
  }
}
