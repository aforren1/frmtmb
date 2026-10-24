# The relative gap between frmtmb's and mclust's log-likelihood on the
# vignette's Old Faithful fit, against the vignette's 1e-7 bound.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.latent))
dd <- data.frame(row = seq_len(nrow(datasets::faithful)))
dd$Y <- data.matrix(datasets::faithful)
f <- frm(bf(Y ~ 1), family = mixture_mvn(K = 2, D = 2), data = dd)
mclustBIC <- mclust::mclustBIC
mc <- mclust::Mclust(dd$Y, G = 2, modelNames = "VVV", verbose = FALSE,
                     control = mclust::emControl(tol = 1e-10))
cat(sprintf("frmtmb %.10f mclust %.10f relative gap %.3e\n",
            as.numeric(logLik(f)), mc$loglik,
            abs(as.numeric(logLik(f)) - mc$loglik) / abs(mc$loglik)))
