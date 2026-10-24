b <- readRDS("dev/review-spatial-log/compat-before.rds")$tab
f <- readRDS("dev/review-spatial-log/compat-after.rds")$tab
for (p in c("car x gr_cov", "car x rr", "car x frm_sample", "gp x spde", "smooth x spde", "equalto x spde")) {
  cat("==", p, "\n  status:", b$status[b$pair == p][1], "\n  BEFORE:", strtrim(b$note[b$pair == p][1], 400),
      "\n  AFTER :", strtrim(f$note[f$pair == p][1], 400), "\n")
}
