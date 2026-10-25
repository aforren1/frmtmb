# predict(newdata = ) on fits whose draw is structured: the refusal's
# words, by model, before and after.
#   Rscript dev/smallfix-predict-newdata.R <lib;lib>
# Log: dev/smallfix-out/predict-newdata-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
say <- function(lab, expr) {
  r <- tryCatch({ force(expr); "answered" }, error = function(e) {
    paste0("[", paste(intersect(class(e), c("frmtmb.latent_error",
                                            "frmtmb_error")), collapse = ","),
           "] ", gsub("\\s+", " ", conditionMessage(e)))
  })
  cat(lab, ":", r, "\n")
}
set.seed(42)
n <- 200
cl <- rbinom(n, 1, 0.4)
x <- rnorm(n)
dd <- data.frame(x = x)
dd$Y <- cbind(ifelse(cl == 1, 0, 3) + 0.5 * x + rnorm(n),
              ifelse(cl == 1, 0, 4) + rnorm(n))
fm <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd)
say("mixture_mvn", predict(fm, newdata = dd[1:5, ]))
set.seed(11)
dh <- do.call(rbind, lapply(1:10, function(id) {
  s <- integer(20); s[1] <- 1L
  for (t in 2:20) s[t] <- sample.int(2, 1, prob = if (s[t - 1] == 1)
    c(.9, .1) else c(.2, .8))
  data.frame(id = id, t = 1:20, x = rnorm(20), y = rnorm(20, c(0, 3)[s]))
}))
fh <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
          data = dh)
say("hmm", predict(fh, newdata = dh[1:5, ]))
set.seed(3)
da <- data.frame(g = factor(rep(1:10, each = 12)), t = rep(1:12, 10),
                 x = rnorm(120))
da$y <- da$x + rnorm(120)
fa <- frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), family = gaussian(),
          data = da)
say("ar(cov = TRUE)", predict(fa, newdata = da[1:5, ]))
fu <- frm(bf(y ~ x + unstr(time = t, gr = g)), family = gaussian(),
          data = da[da$t <= 3, ])
say("unstr()", predict(fu, newdata = da[1:5, ]))
# the guard absent: an unstructured draw on newdata answers
fg <- frm(bf(y ~ x), family = gaussian(), data = da)
say("gaussian", predict(fg, newdata = da[1:5, ]))
