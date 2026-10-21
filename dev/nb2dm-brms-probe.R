# What brms 2.23.0 does with negbinomial2 and dirichlet_multinomial:
# Stan code, default priors, the family's slots, and which variants and
# spellings exist. No compile.
#   Rscript dev/nb2dm-brms-probe.R > dev/nb2dm-out/brms-probe.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 40
d <- data.frame(x = rnorm(n), g = gl(4, 10))
d$y <- rnbinom(n, mu = exp(1 + 0.3 * d$x), size = 2)
d$z <- rnorm(n)

cat("\n== negbinomial2 exported? ", "negbinomial2" %in% getNamespaceExports("brms"), "\n")
f2 <- brmsfamily("negbinomial2")
print(f2)
str(unclass(f2)[c("family", "link", "dpars", "link_sigma", "type")])
cat("\n== negbinomial2 variants in brms's family list\n")
fl <- grep("negbinomial|multinomial", ls(asNamespace("brms"), all.names = TRUE),
           value = TRUE)
print(grep("^[.]family_", fl, value = TRUE))
cat("\n== stancode y ~ x, negbinomial2\n")
cat(stancode(y ~ x, data = d, family = brmsfamily("negbinomial2")))
cat("\n== default_prior y ~ x, negbinomial2\n")
print(default_prior(y ~ x, data = d, family = brmsfamily("negbinomial2")))
cat("\n== default_prior bf(y ~ x, sigma ~ z), negbinomial2\n")
print(default_prior(bf(y ~ x, sigma ~ z), data = d,
                    family = brmsfamily("negbinomial2")))
cat("\n== default_prior y ~ x, negbinomial (control)\n")
print(default_prior(y ~ x, data = d, family = negbinomial()))
cat("\n== default_prior bf(y ~ x, shape ~ z), negbinomial (control)\n")
print(default_prior(bf(y ~ x, shape ~ z), data = d, family = negbinomial()))
cat("\n== stancode bf(y ~ x, sigma ~ z), negbinomial2\n")
cat(stancode(bf(y ~ x, sigma ~ z), data = d,
             family = brmsfamily("negbinomial2")))
cat("\n== stancode y | rate(z2) ~ x, negbinomial2\n")
d$z2 <- exp(d$z)
cat(stancode(y | rate(z2) ~ x, data = d, family = brmsfamily("negbinomial2")))
cat("\n== stancode y | cens(c) ~ x, negbinomial2\n")
d$c <- rep(c(0, 1), length.out = n)
cat(stancode(y | cens(c) ~ x, data = d, family = brmsfamily("negbinomial2")))
cat("\n== stancode y | trunc(lb = 1) ~ x, negbinomial2\n")
cat(stancode(y | trunc(lb = 1) ~ x, data = d[d$y >= 1, ],
             family = brmsfamily("negbinomial2")))
cat("\n== negbinomial2(sqrt), link_sigma identity\n")
print(tryCatch(brmsfamily("negbinomial2", link = "sqrt",
                          link_sigma = "identity"),
               error = function(e) conditionMessage(e)))
print(tryCatch(brmsfamily("negbinomial2", link_sigma = "softplus"),
               error = function(e) conditionMessage(e)))
cat("\n== mixture(negbinomial2, poisson)\n")
print(tryCatch(stancode(y ~ x, data = d,
                        family = mixture(brmsfamily("negbinomial2"),
                                         poisson())),
               error = function(e) conditionMessage(e)))

# dirichlet_multinomial
Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, c(0.2, 0.3, 0.5))))
colnames(Y) <- c("a", "b", "c")
d$Y <- Y
d$tr <- rowSums(Y)
cat("\n== dirichlet_multinomial()\n")
fd <- dirichlet_multinomial()
print(fd)
str(unclass(fd)[c("family", "link", "dpars", "link_phi", "refcat", "type")])
cat("\n== stancode Y | trials(tr) ~ x, dirichlet_multinomial\n")
cat(stancode(Y | trials(tr) ~ x, data = d, family = dirichlet_multinomial()))
cat("\n== default_prior Y | trials(tr) ~ x, dirichlet_multinomial\n")
print(default_prior(Y | trials(tr) ~ x, data = d,
                    family = dirichlet_multinomial()))
cat("\n== default_prior bf(Y | trials(tr) ~ x, phi ~ z), dm\n")
print(default_prior(bf(Y | trials(tr) ~ x, phi ~ z), data = d,
                    family = dirichlet_multinomial()))
cat("\n== stancode without trials()\n")
print(tryCatch(stancode(Y ~ x, data = d, family = dirichlet_multinomial()),
               error = function(e) conditionMessage(e)))
cat("\n== standata trials vs rowSums mismatch\n")
d2 <- d
d2$tr[1] <- d2$tr[1] + 1
print(tryCatch(names(standata(Y | trials(tr) ~ x, data = d2,
                              family = dirichlet_multinomial())),
               error = function(e) conditionMessage(e)))
cat("\n== refcat = 'b', refcat = NA\n")
print(tryCatch(default_prior(Y | trials(tr) ~ x, data = d,
                             family = dirichlet_multinomial(refcat = "b")),
               error = function(e) conditionMessage(e)))
print(tryCatch(default_prior(Y | trials(tr) ~ x, data = d,
                             family = dirichlet_multinomial(refcat = NA)),
               error = function(e) conditionMessage(e)))
cat("\n== mvbf with dirichlet_multinomial and gaussian\n")
print(tryCatch({
  sc <- stancode(mvbf(bf(Y | trials(tr) ~ x, family = dirichlet_multinomial()),
                      bf(z ~ x, family = gaussian())), data = d)
  cat(sc)
  "ok"
}, error = function(e) conditionMessage(e)))
cat("\n== mvbf with multinomial and gaussian (control)\n")
print(tryCatch({
  sc <- stancode(mvbf(bf(Y | trials(tr) ~ x, family = multinomial()),
                      bf(z ~ x, family = gaussian())), data = d)
  "ok"
}, error = function(e) conditionMessage(e)))
cat("\n== random intercept dm\n")
cat(stancode(Y | trials(tr) ~ x + (1 | g), data = d,
             family = dirichlet_multinomial()))
cat("\n== mixture(dirichlet_multinomial)\n")
print(tryCatch(stancode(Y | trials(tr) ~ x, data = d,
                        family = mixture(dirichlet_multinomial(),
                                         dirichlet_multinomial())),
               error = function(e) conditionMessage(e)))
cat("\n== dm link other than logit\n")
print(tryCatch(dirichlet_multinomial(link = "probit"),
               error = function(e) conditionMessage(e)))
print(tryCatch(dirichlet_multinomial(link_phi = "identity"),
               error = function(e) conditionMessage(e)))
print(tryCatch(dirichlet_multinomial(link_phi = "softplus"),
               error = function(e) conditionMessage(e)))
cat("\n== negbinomial links accepted for sigma\n")
print(brms:::.family_negbinomial2()$links)
print(brms:::valid_dpars(brmsfamily("negbinomial2")))
print(brms:::valid_dpars(dirichlet_multinomial()))
