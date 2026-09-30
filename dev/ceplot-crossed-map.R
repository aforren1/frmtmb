# Lane ceplot: how brms 2.23.0 numbers the three group terms of
# y ~ x + (1 | g) + (1 | h) + (1 | g:h) and orders their levels, so that
# dev/ceplot-crossed-brms.R can initialize brms at frmtmb's draws.
# Data seed 49, as dev/postfit2-p2-checks.R's P1-M2.
#   Rscript dev/ceplot-crossed-map.R > dev/ceplot-log/crossed-map.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
fc <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = dc)
cat("frmtmb labels:\n")
print(frmtmb::brms_par_labels(fc))
sd <- brms::standata(y ~ x + (1 | g) + (1 | h) + (1 | g:h), data = dc)
for (k in 1:3) {
  cat("N_", k, " = ", sd[[paste0("N_", k)]], ", M_", k, " = ",
      sd[[paste0("M_", k)]], "\n", sep = "")
}
init <- list(list(b = array(0.5, 1), Intercept = 1, sigma = 1,
                  sd_1 = array(1, 1), sd_2 = array(2, 1), sd_3 = array(3, 1),
                  z_1 = matrix(seq_len(sd$N_1) / 100, 1),
                  z_2 = matrix(seq_len(sd$N_2) / 100, 1),
                  z_3 = matrix(seq_len(sd$N_3) / 100, 1)))
b <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g) + (1 | h) + (1 | g:h), data = dc,
      algorithm = "fixed_param", chains = 1, iter = 1, warmup = 0,
      init = init, refresh = 0, seed = 1, silent = 2)))
m <- as_draws_matrix(b)
print(round(m[1, grepl("^sd_|^r_", colnames(m))], 4))
cat("done\n")
