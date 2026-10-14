# Lane optstart: what a frame holds for the start code (linpreds, aterms)
#   Rscript dev/optstart-peek.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(50), g = factor(rep(1:10, 5)), w = runif(50),
                o = rnorm(50))
d$y <- 1 + 2 * d$x + rnorm(10)[d$g] + rnorm(50)
f <- frm(bf(y | weights(w) ~ x + offset(o) + (1 | g), sigma ~ x),
         data = d, family = gaussian())
fr <- f$frame
cat("frame names:", names(fr), "\n")
for (k in names(fr$linpreds)) {
  lp <- fr$linpreds[[k]]
  cat("linpred", k, ":", names(lp), "\n")
  cat("   par", lp$par, "dpar", lp$dpar, "resp", lp$resp, "idx",
      lp$idx, "link", lp$link$name, "\n")
}
str(fr$aterm_values, max.level = 2)
str(fr$par_template)
cat("re_blocks:", length(fr$re_blocks), "\n")
str(fr$re_blocks[[1]][c("dim", "covstruct", "theta_idx", "components")],
    max.level = 3)
print(f$start %||% "no start slot")
print(names(f))
