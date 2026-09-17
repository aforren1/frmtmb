# The nonlinear body of brms's brmsfit_example2 fails to evaluate in frm().
# Which piece of the body does it? Rscript dev/brmsport-probe-nl.R
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib", "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5)
d <- data.frame(Trt = rep(0:1, 20), Age = rnorm(40))
d$count <- rgamma(40, shape = 5, rate = 5 / (1 / (1 + exp(-(0.5 + 0.2 * d$Age))) * exp(0.3 * d$Trt)))
try1 <- function(lbl, f, fam = Gamma("identity"), ...) {
  r <- tryCatch({fit <- suppressWarnings(frm(f, d, family = fam, ...)); sprintf("OK logLik %.4f", logLik(fit))},
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-40s %s\n", lbl, r))
}
try1("1/(1 + exp(-a)) * exp(b * Trt)", bf(count ~ 1/(1 + exp(-a)) * exp(b * Trt), a ~ Age, b ~ 1, nl = TRUE))
try1("plogis-free: exp(b*Trt)/(1+exp(-a))", bf(count ~ exp(b * Trt) / (1 + exp(-a)), a ~ Age, b ~ 1, nl = TRUE))
try1("1/(1 + exp(a))", bf(count ~ 1/(1 + exp(a)), a ~ Age, nl = TRUE))
try1("1/(1 + exp(0 - a))", bf(count ~ 1/(1 + exp(0 - a)), a ~ Age, nl = TRUE))
try1("exp(-a)", bf(count ~ exp(-a), a ~ Age, nl = TRUE))
try1("-a + 30", bf(count ~ -a + 30, a ~ Age, nl = TRUE), fam = gaussian())
try1("1/a", bf(count ~ 1/a, a ~ Age, nl = TRUE), fam = gaussian())
try1("a * 1", bf(count ~ a * 1, a ~ Age, nl = TRUE), fam = gaussian())
cat("-- with brms's example data\n")
d <- as.data.frame(get("brmsfit_example2", asNamespace("brms"))$data)
try1("ex2 as is", bf(count | weights(AgeSD) ~ 1/(1 + exp(-a)) * exp(b * Trt), a ~ Age + (1 | ID1 | patient), b ~ Age + (1 | ID1 | patient), nl = TRUE))
try1("ex2 no weights", bf(count ~ 1/(1 + exp(-a)) * exp(b * Trt), a ~ Age + (1 | ID1 | patient), b ~ Age + (1 | ID1 | patient), nl = TRUE))
try1("ex2 no RE", bf(count ~ 1/(1 + exp(-a)) * exp(b * Trt), a ~ Age, b ~ Age, nl = TRUE))
d$Trtn <- as.numeric(as.character(d$Trt))
try1("ex2 numeric Trt", bf(count | weights(AgeSD) ~ 1/(1 + exp(-a)) * exp(b * Trtn), a ~ Age + (1 | ID1 | patient), b ~ Age + (1 | ID1 | patient), nl = TRUE))
try1("factor Trt, gaussian, a*Trt", bf(count ~ a * Trt, a ~ 1, nl = TRUE), fam = gaussian())
d2 <- d; d2$Trt <- d2$Trtn
d2$Trt <- factor(d2$Trt)
cat("brms, same data: C matrix of the factor\n")
sd <- brms::make_standata(brms::bf(count | weights(AgeSD) ~ 1/(1 + exp(-a)) * exp(b * Trt), a ~ Age + (1 | ID1 | patient), b ~ Age + (1 | ID1 | patient), nl = TRUE), data = d, family = brms::brmsfamily("gamma", "identity"))
print(head(cbind(sd$C_1, Trt = d$Trt)))
