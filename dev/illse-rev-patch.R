f <- "C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-glm.R"
x <- readLines(f)
i <- grep("out$marginaleffects <- mf$v$std.error", x, fixed = TRUE)
x[i] <- paste0(
"  if (!inherits(mf$v, \"rev_err\")) {\n",
"    est_f <- if (fam == \"binom\") fr$v[, \"Estimate\"] else as.vector(fr$v[, \"Estimate\", ])\n",
"    se_f <- if (fam == \"binom\") fr$v[, \"Estimate\"] * 0 + sqrt(vr) else as.vector(sqrt(vr))\n",
"    est_m <- mf$v$estimate\n",
"    j <- vapply(est_m, function(v) which.min(abs(est_f - v)), 1L)\n",
"    cat(\"  mfx est match rel\", sprintf(\"%.1e\", max(abs(est_m / est_f[j] - 1))), \"\n\")\n",
"    out$marginaleffects <- mf$v$std.error / se_f[j]\n",
"  }")
writeLines(x, f)
