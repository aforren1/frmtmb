# Lane bigx: the grid of the brief. Offsets m = 1e3 to 1e8 x families
# (poisson, bernoulli, gaussian, Gamma with log link) x with and without
# (1 | g) x data seeds. Each cell fits y ~ w (+ (1 | g)) with w = z + m
# and compares it with the same model on z, an exact reparametrization:
# the slope and the log-likelihood must agree. lme4 2.1-0 (or glm())
# fits the same w as a second reference.
#   BIGX_ARM=before|after Rscript dev/bigx-grid.R <out.tsv> [seeds]
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
suppressPackageStartupMessages(library(lme4))
cat("lme4", format(packageVersion("lme4")), "\n")
a <- commandArgs(TRUE)
out <- a[1]
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:5
offs <- c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
fams <- c("poisson", "binomial", "gaussian", "Gamma")
rows <- list()
lme4_slope <- function(fam, d, re) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch({
    if (re) {
      f <- switch(fam,
        gaussian = lmer(y ~ w + (1 | g), data = d, REML = FALSE),
        poisson = glmer(y ~ w + (1 | g), data = d, family = poisson()),
        binomial = glmer(y ~ w + (1 | g), data = d, family = binomial()),
        Gamma = glmer(y ~ w + (1 | g), data = d,
                      family = Gamma(link = "log")))
      fe <- fixef(f)
      c(slope = unname(fe["w"]), ll = as.numeric(logLik(f)))
    } else {
      f <- glm(y ~ w, data = d, family = switch(fam,
        gaussian = gaussian(), poisson = poisson(), binomial = binomial(),
        Gamma = Gamma(link = "log")))
      c(slope = unname(coef(f)["w"]), ll = NA_real_)
    }
  }, error = function(e) c(slope = NA_real_, ll = NA_real_)),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleMessage")
  })
  list(slope = r[["slope"]], ll = r[["ll"]], nw = length(w))
}
t_all <- proc.time()[["elapsed"]]
for (fam in fams) for (re in c(TRUE, FALSE)) for (m in offs) for (s in seeds) {
  d <- bigx_data(fam, m, s, re = re)
  fm <- bigx_family(fam)
  rhs <- if (re) " + (1 | g)" else ""
  fc <- bigx_fit(bf(as.formula(paste0("y ~ z", rhs))), family = fm, data = d)
  t0 <- proc.time()[["elapsed"]]
  fw <- bigx_fit(bf(as.formula(paste0("y ~ w", rhs))), family = fm, data = d)
  tw <- proc.time()[["elapsed"]] - t0
  lm4 <- lme4_slope(fam, d, re)
  sec <- if (!fc$err) bigx_se(fc$f) else NA
  row <- data.frame(family = fam, re = re, m = m, seed = s,
                    c_err = fc$err, c_warn = length(fc$w),
                    c_slope = if (fc$err) NA else fc$f$opt$par[[2]],
                    c_se = sec[2],
                    w_err = fw$err, w_code = NA, w_warn = length(fw$w),
                    w_classes = paste(unique(fw$wc), collapse = ","),
                    w_first = if (length(fw$w)) gsub("\t", " ", fw$w[1]) else "",
                    w_haslope = NA, slope_z = NA, se_ratio = NA, dll = NA,
                    theta_z = NA, w_sec = tw,
                    lme4_slope_z = (lm4$slope - (if (fc$err) NA else
                                                   fc$f$opt$par[[2]])) / sec[2],
                    lme4_dll = NA, lme4_nw = lm4$nw)
  if (fw$err) {
    row$w_first <- substr(gsub("[\t\n]", " ", conditionMessage(fw$f)), 1, 90)
  } else if (!fc$err) {
    f <- fw$f
    b <- f$opt$par
    row$w_code <- f$opt$convergence
    row$w_haslope <- length(b) == length(fc$f$opt$par)
    llw <- as.numeric(logLik(f))
    row$dll <- llw - as.numeric(logLik(fc$f))
    if (isTRUE(row$w_haslope)) {
      se <- bigx_se(f)
      row$slope_z <- (b[[2]] - fc$f$opt$par[[2]]) / sec[2]
      row$se_ratio <- se[2] / sec[2]
      if (length(b) >= 3) {
        row$theta_z <- (b[[length(b)]] - fc$f$opt$par[[length(b)]]) /
          sec[length(b)]
      }
    }
    row$lme4_dll <- lm4$ll - as.numeric(logLik(fc$f))
  }
  rows[[length(rows) + 1L]] <- row
  cat(sprintf("%-8s re %-5s m %-6g s %d | slope z %+9.3f SE %.4f dll %+.2e code %s warn %d %s | lme4 z %+.3f\n",
              fam, re, m, s, row$slope_z, row$se_ratio, row$dll,
              format(row$w_code), row$w_warn, row$w_classes,
              row$lme4_slope_z))
}
tab <- do.call(rbind, rows)
write.table(tab, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("total seconds", proc.time()[["elapsed"]] - t_all, "\n")
