# Short check 2, items 1 and 3: summaries of review-latres-tk.R.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/tk"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE, comment.char = "")))
x <- x[is.finite(x$ll), ]
x$rs <- exp(x$ls - x$x_ls)
x$dar <- x$ar - x$x_ar
x$cov <- abs(x$ls - log(x$s_true)) <= 1.96 * x$se_ls
# one row per design: bad by the lane's definition
des <- do.call(rbind, lapply(split(x, paste(x$set, x$design)), function(d) {
  data.frame(set = d$set[1], design = d$design[1], n = nrow(d), ratio = median(d$rs, na.rm = TRUE), dar = median(d$dar, na.rm = TRUE),
             cov = mean(d$cov, na.rm = TRUE), warned = sum(d$warned == 1, na.rm = TRUE), refused = sum(d$margin < 2, na.rm = TRUE),
             code1 = sum(d$code != 0, na.rm = TRUE), boxed = sum(d$boxed == 1, na.rm = TRUE),
             Cn = median(abs(d$C) / d$n, na.rm = TRUE), bls = median(abs(d$b_ls), na.rm = TRUE), bar = median(abs(d$b_ar), na.rm = TRUE),
             Cgap = median((abs(d$C_up - d$C) + abs(d$C_dn - d$C)), na.rm = TRUE),
             Cci = median(abs(d$C_chi - d$C_clo), na.rm = TRUE),
             ttk = mean(d$t_tk, na.rm = TRUE))
}))
des$bad <- with(des, ratio < .9 | ratio > 1.15 | abs(dar) > .15 | cov < .8)
x <- merge(x, des[, c("set", "design", "bad")], by = c("set", "design"))
x$st <- with(x, pmax(abs(C_chi - C), abs(C_clo - C), na.rm = TRUE))
x$tk_b <- with(x, pmax(abs(b_ls), abs(b_ar), na.rm = TRUE))
x$Cn <- abs(x$C) / x$n
x$Cd2 <- with(x, pmax(abs(C_up - C), abs(C_dn - C), na.rm = TRUE))
cat("== the current warning, by set: hit among bad designs' fits, false alarm among others\n")
rate <- function(fire, lab) {
  o <- do.call(rbind, lapply(split(seq_len(nrow(x)), x$set), function(i) {
    b <- x$bad[i]; w <- fire[i]
    data.frame(set = x$set[i][1], hit = round(mean(w[b], na.rm = TRUE), 3), fa = round(mean(w[!b], na.rm = TRUE), 3), nb = sum(b), ng = sum(!b))
  }))
  all <- data.frame(set = "ALL", hit = round(mean(fire[x$bad], na.rm = TRUE), 3), fa = round(mean(fire[!x$bad], na.rm = TRUE), 3), nb = sum(x$bad), ng = sum(!x$bad))
  o <- rbind(o, all)
  cat(sprintf("%-34s %s\n", lab, paste(sprintf("%s %.3f/%.3f", o$set, o$hit, o$fa), collapse = "  ")))
}
rate(x$warned == 1, "lane's condition (hit/false alarm)")
rate((x$warned == 1) | (x$margin < 2 & !is.na(x$margin)), "lane's condition + margin refusal")
for (k in c(.0005, .001, .002, .005, .01)) rate(x$Cn > k, sprintf("|C|/n > %.4f", k))
for (k in c(.05, .1, .2, .3)) rate(x$tk_b > k | is.na(x$tk_b), sprintf("|shift| > %.2f or NA", k))
for (k in c(.2, .3, .5, 1, 1.5, 2)) rate(x$Cd2 > k | is.na(x$Cd2), sprintf("C moves > %.1f at sderr x2 or /2", k))
for (k in c(.1, .2, .3, .5, 1)) rate(x$st > k | is.na(x$st), sprintf("C moves > %.1f across the CI", k))
cat("\n== each design\n")
print(des[order(des$set, des$ratio), c("set", "design", "n", "ratio", "dar", "cov", "bad", "warned", "refused", "code1", "boxed", "Cn", "bls", "bar", "Cgap", "Cci")], digits = 2, row.names = FALSE)
cat("\nmean seconds of the correction with its gradient, all fits:", round(mean(x$t_tk, na.rm = TRUE), 3), " largest n", max(x$n, na.rm = TRUE), "\n")
cat("\n== does Laplace + C read the exact value at the fit's point?\n")
y <- x[is.finite(x$x_at) & is.finite(x$C) & abs(x$x_at) < 1e6, ]
y$g0 <- y$ll - y$x_at
y$g1 <- y$ll + y$C - y$x_at
print(do.call(rbind, lapply(split(y, y$set), function(d) data.frame(set = d$set[1], fits = nrow(d), med_abs_laplace = round(median(abs(d$g0)), 3), med_abs_corrected = round(median(abs(d$g1)), 3), q90_laplace = round(quantile(abs(d$g0), .9), 3), q90_corrected = round(quantile(abs(d$g1), .9), 3), row.names = NULL))))
cat("\n== does the shift move sderr toward the exact MLE?\n")
z <- x[is.finite(x$b_ls) & abs(x$ls) < 10 & abs(x$b_ls) < 3, ]
print(do.call(rbind, lapply(split(z, z$set), function(d) data.frame(set = d$set[1], fits = nrow(d), med_err_before = round(median(abs(d$ls - d$x_ls)), 3), med_err_after = round(median(abs(d$ls + d$b_ls - d$x_ls)), 3)))))
cat("\n== item 1: the benign link-scale designs on the final build\n")
link <- c("ginv_b1","ginv_b2","ginv_b3","ginv_s5","gident_b","gid_s5","ig_b1","ig_b2","ig_s5","pid_b","pid_20","pid_8","bid_b","gid_l","bid_l","ginv_l","pid_l")
L <- x[x$set == "m4b" & x$design %in% link, ]
L$ok <- abs(L$ls - L$x_ls) < log(1.25) & abs(L$ar - L$x_ar) < .3 & abs(L$ll - L$x_at) < 1
L$ref <- L$margin < 2
cat("fits", nrow(L), " refused", sum(L$ref), " of them read true", sum(L$ok[L$ref]), " kept", sum(!L$ref), " of them read true", sum(L$ok[!L$ref]), round(mean(L$ok[!L$ref]), 3), " warned among kept", sum(L$warned[!L$ref] == 1), " boxed", sum(L$boxed == 1), "\n")
print(do.call(rbind, lapply(split(L, L$design), function(d) data.frame(design = d$design[1], refused = sum(d$ref), kept = sum(!d$ref), kept_true = sum(d$ok[!d$ref]), kept_warned = sum(d$warned[!d$ref] == 1), med_margin = round(median(d$margin), 2)))), row.names = FALSE)
