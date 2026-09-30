# hypothesis() on the reviewer's equated hmm (seed 4001), naming the
# equated cell and its target
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
set.seed(4001)
G <- 20; Tn <- 25
z <- unlist(lapply(seq_len(G), function(g) {
  s <- integer(Tn); s[1] <- sample(1:2, 1)
  for (t in 2:Tn) s[t] <- if (runif(1) < 0.8) s[t - 1] else 3 - s[t - 1]
  s
}))
dd <- data.frame(id = rep(seq_len(G), each = Tn), t = rep(seq_len(Tn), G))
dd$y <- rnorm(nrow(dd), c(0, 3)[z], 0.7)
h1 <- suppressWarnings(frm(bf(y ~ 1, tr12 = "tr22"),
          family = hmm(K = 2, gaussian(), time = t, group = id,
                       init = "stationary"), data = dd))
print(tryCatch(hypothesis(h1, "tr12_Intercept = 0"),
               error = function(e) conditionMessage(e)))
print(hypothesis(h1, "tr22_Intercept = 0")$hypothesis)
print(summary(h1)$spec_pars)
