# Reviewer: what brms 2.23.0 builds for the parser cases of
# dev/famsimplex-rev-parse.R. Output: dev/famsimplex-rev-brmsparse.log
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:15, length.out = n)))
E <- cbind(0, 0.3 + 0.5 * d$z, -0.2 + 0.6 * d$x)
P <- exp(E) / rowSums(exp(E))
d$f <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))])
d$N <- 10L
d$Ym <- t(vapply(seq_len(n), function(i) rmultinom(1, 10, P[i, ])[, 1],
                 numeric(3)))
colnames(d$Ym) <- c("a", "b", "c")
d$y <- rnorm(n)
sc <- function(lab, f, fam, dat = d) {
  cat("\n---", lab, "\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch({
    s <- standata(f, data = dat, family = fam)
    code <- stancode(f, data = dat, family = fam)
    xs <- grep("^X", names(s), value = TRUE)
    for (x in xs) cat(x, ":", colnames(s[[x]]), "\n")
    ln <- grep("mu[0-9a-z]* (\\+|=)", strsplit(code, "\n")[[1]], value = TRUE)
    cat(ln, sep = "\n")
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e))),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  if (length(w)) cat("WARN:", unique(w), sep = "\n  ")
  cat(r, "\n")
}
mix <- mixture(gaussian, gaussian)
sc("mixture y ~ x, mu1 ~ z", bf(y ~ x, mu1 ~ z), mix)
sc("mixture y ~ w, mu1 ~ x, mu2 ~ z", bf(y ~ w, mu1 ~ x, mu2 ~ z), mix)
sc("mixture y ~ x, mu1 = 0", bf(y ~ x, mu1 = 0), mix)
sc("categorical f ~ w, mub ~ z, muc ~ x", bf(f ~ w, mub ~ z, muc ~ x),
   categorical())
sc("categorical f ~ 1, mub ~ z, muc ~ x", bf(f ~ 1, mub ~ z, muc ~ x),
   categorical())
sc("multinomial mub ~ z", bf(Ym | trials(N) ~ x, mub ~ z), multinomial())
sc("categorical f ~ x, mu ~ z", bf(f ~ x, mu ~ z), categorical())
d2 <- d
d2$w[1:10] <- NA
cat("\n--- NA in w only, f ~ w, mub ~ z, muc ~ x: N =",
    tryCatch(standata(bf(f ~ w, mub ~ z, muc ~ x), data = d2,
                      family = categorical())$N,
             error = function(e) conditionMessage(e)), "\n")
