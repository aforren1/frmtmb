# Summary of review-latres-m1b.R
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m1b"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
tr <- c(zip=.5,zip_hi=.5,zip_lo=.5,hp=.5,zib=.6,cum3=1,cum5=1,cum7=1,nb_m1=.5,nb_m5=.5,nb_m20=.5,nb_s1=1,beta_lo=.5,binom4=.8,ln_s2=1,ln_s05=.5,sn_u1=.5,sn_u10=5)
b1t <- c(zip=.4,zip_hi=.4,zip_lo=.4,hp=.4,zib=.4,cum3=.5,cum5=.5,cum7=.5,nb_m1=.4,nb_m5=.4,nb_m20=.4,nb_s1=.4,beta_lo=.4,binom4=.5,ln_s2=.4,ln_s05=.4,sn_u1=.4,sn_u10=4)
out <- do.call(rbind, lapply(split(x, x$design), function(d) {
  dn <- d$design[1]
  rs <- exp(d$ls) / exp(d$x_ls)
  data.frame(design = dn, n = nrow(d), errs = sum(nzchar(d$err) & !is.na(d$err)), code1 = sum(d$code != 0, na.rm = TRUE),
    info = median(d$info, na.rm = TRUE), warned = sum(d$warned == 1, na.rm = TRUE),
    s_true = tr[dn], s_frm = median(exp(d$ls), na.rm = TRUE), s_ex = median(exp(d$x_ls), na.rm = TRUE), ratio = median(rs, na.rm = TRUE),
    ar_frm = median(d$ar, na.rm = TRUE), ar_ex = median(d$x_ar, na.rm = TRUE),
    b1_frm = median(d$b1, na.rm = TRUE), b1_ex = median(d$x_b1, na.rm = TRUE),
    cov_ls = mean(abs(d$ls - log(tr[dn])) <= 1.96 * d$se_ls, na.rm = TRUE),
    cov_b1 = mean(abs(d$b1 - b1t[dn]) <= 1.96 * d$se_b1, na.rm = TRUE),
    bad = abs(median(d$ar - d$x_ar, na.rm = TRUE)) > .15 | median(rs, na.rm = TRUE) < .9 | median(rs, na.rm = TRUE) > 1.15, row.names = NULL)
}))
print(out[order(out$info), ], digits = 3, row.names = FALSE)
