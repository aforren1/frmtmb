# Punch round 2: the round-1 rule (reach > 50 or collapsed) on the
# margin designs' degenerate fits, probit saturation apart, for the
# findings table (from dev/ordmix-p2-crit-log/).
dir <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-p2-crit-log"
n <- c(probit = 0, other = 0)
w <- n
for (f in list.files(dir, "^margin_.*[.]txt$", full.names = TRUE)) {
  for (l in grep("label=degenerate", readLines(f), value = TRUE)) {
    kv <- strsplit(l, " ")[[1]]
    v <- setNames(sub("^[^=]*=", "", kv), sub("=.*$", "", kv))
    r <- max(as.numeric(strsplit(v[["reach"]], ",")[[1]]))
    col <- any(strsplit(v[["collapsed"]], ",")[[1]] == "TRUE")
    k <- if (grepl("probit", v[["cfg"]])) "probit" else "other"
    n[k] <- n[k] + 1
    w[k] <- w[k] + (r > 50 || col)
  }
}
print(rbind(fits = n, round1_warns = w))
