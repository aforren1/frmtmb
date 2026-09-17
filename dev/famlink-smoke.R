.libPaths(c("C:/Users/adf44/source/r/famlink-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(find.package("frmtmb"), "\n")
print(weibull())
print(cumulative("probit"))
str(unclass(beta_binomial()), max.level = 1, give.attr = FALSE)
print(try(student()$lin))
print(student()$nope)
print(try({f <- student(); f$link <- "log"}))
f <- student()
f[["links"]][["sigma"]] <- "softplus"
cat("after [[<-: link_sigma =", f$link_sigma, "\n")
f$links$nu <- "identity"
cat("after $<-: link_nu =", f$link_nu, "\n")
m <- mixture(gaussian(), student())
print(names(unclass(m)))
cat("mixture link:", m$link, "\n")
k <- categorical(levels = c("a", "b", "c"))
print(names(unclass(k))); cat("categorical link:", k$link, "\n")
mn <- multinomial(3); print(grep("^link", names(unclass(mn)), value = TRUE))
print(brmsfamily("normal")$family)
print(brmsfamily("Gamma", identity)$link)
print(try(brmsfamily(c("poisson", "log"))))
print(try(frm_family(c("poisson", "log"), link = "sqrt")))
bf <- brms::brmsfamily("beta_binomial", link = "probit", link_phi = "softplus")
fb <- frmtmb:::as_frmtmb_family(bf)
cat("brms family through as_frmtmb_family: link", fb$link, "link_phi", fb$link_phi, "\n")
set.seed(3)
d <- data.frame(y = rpois(40, 3), x = rnorm(40))
print(fixef(frm(y ~ x, d, family = "zi_poisson"))$mu)
print(try(mixture(gaussian, hurdle_gamma)))
if (requireNamespace("frmtmb.eam", quietly = TRUE)) {
  w <- frmtmb.eam::wiener()
  print(grep("^link", names(unclass(w)), value = TRUE))
  print(unlist(unclass(w)[grep("^link", names(unclass(w)))][-c(2, 3)]))
}
if (requireNamespace("frmtmb.learn", quietly = TRUE)) {
  r <- try(frmtmb.learn::rlddm)
}
