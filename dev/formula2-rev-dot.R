# Reviewer: try to falsify the `.` expansion claims against brms 2.23.0's
# validate_formula()/standata(), and probe update() of a dot fit.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
show <- function(x) {
  if (inherits(x, "ERR")) paste("ERR:", substr(x, 1, 150)) else x
}
set.seed(6)
n <- 120
d <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n),
                g = factor(rep(1:6, each = 20)))
d$y <- d$y + d$x1

xcols <- function(fr, key) {
  cn <- colnames(fr$linpreds[[key]]$X)
  cn[cn == "(Intercept)"] <- "Intercept"
  cn
}
cmp <- function(label, ff, fb = ff, data = d, keys = "y.mu",
                sdn = "X", fam = NULL, bfam = NULL) {
  fr <- q(frm(ff, data = data, family = fam, dry_run = "frame"))
  sd <- q(if (is.null(bfam)) brms::standata(fb, data) else
    brms::standata(fb, data, family = bfam))
  cat(sprintf("%-30s\n", label))
  for (i in seq_along(keys)) {
    a <- if (inherits(fr, "ERR")) show(fr) else xcols(fr, keys[i])
    b <- if (inherits(sd, "ERR")) show(sd) else colnames(sd[[sdn[i]]])
    cat("   ", keys[i], "frm [", paste(a, collapse = ","), "] brms [",
        paste(b, collapse = ","), "]", identical(a, b), "\n")
  }
  invisible(fr)
}
B <- brms::bf
cmp("y ~ .", y ~ .)
cmp("y ~ . - x1", y ~ . - x1)
cmp("y ~ .^2", y ~ .^2)
cmp("y ~ . + (1 | g)", y ~ . + (1 | g))
cmp("y ~ . - g + (1 | g)", y ~ . - g + (1 | g))
cmp("y ~ offset(x2) + .", y ~ offset(x2) + .)
ds <- d; ds$sigma <- rnorm(n)
cmp("data column named sigma", y ~ ., data = ds)
cmp("sigma ~ .", bf(y ~ x1, sigma ~ .), B(y ~ x1, sigma ~ .),
    keys = c("y.mu", "y.sigma"), sdn = c("X", "X_sigma"))
cmp("nl a ~ .", bf(y ~ a + b * x1, a ~ ., b ~ 1, nl = TRUE),
    B(y ~ a + b * x1, a ~ ., b ~ 1, nl = TRUE), keys = "y.a", sdn = "X_a")

cat("== multivariate\n")
dm <- d; dm$y2 <- rnorm(n)
cmp("bf(y ~ .) + bf(y2 ~ .)", bf(y ~ .) + bf(y2 ~ .) + set_rescor(FALSE),
    B(y ~ .) + B(y2 ~ .) + brms::set_rescor(FALSE), data = dm,
    keys = c("y.mu", "y2.mu"), sdn = c("X_y", "X_y2"))
cmp("mvbind(y, y2) ~ .", bf(mvbind(y, y2) ~ .) + set_rescor(FALSE),
    B(brms::mvbind(y, y2) ~ .) + brms::set_rescor(FALSE), data = dm,
    keys = c("y.mu", "y2.mu"), sdn = c("X_y", "X_y2"))
cmp("mvbind(y, y2) ~ . (rescor)", bf(mvbind(y, y2) ~ .),
    B(brms::mvbind(y, y2) ~ .), data = dm,
    keys = c("y.mu", "y2.mu"), sdn = c("X_y", "X_y2"))

cat("== addition terms\n")
dw <- d; dw$w <- runif(n, 0.5, 1.5); dw$cc <- rbinom(n, 1, 0.1)
cmp("y | weights(w) ~ .", y | weights(w) ~ ., data = dw)
cmp("y | cens(cc) ~ .", y | cens(cc) ~ ., data = dw)
dt <- d; dt$k <- rbinom(n, 10, 0.5); dt$nt <- 10L
cmp("k | trials(nt) ~ .", bf(k | trials(nt) ~ . - y), data = dt,
    fam = binomial(), bfam = brms::brmsfamily("binomial"), keys = "k.mu")

cat("== entry points\n")
cat("get_prior y ~ . identical to written out:",
    identical(q(get_prior(y ~ ., data = d)),
              q(get_prior(y ~ x1 + x2 + g, data = d))), "\n")
cat("default_prior y ~ . identical to written out:",
    identical(q(default_prior(y ~ ., data = d)),
              q(default_prior(y ~ x1 + x2 + g, data = d))), "\n")
cat("par_template identical:",
    identical(q(par_template(y ~ ., data = d)),
              q(par_template(y ~ x1 + x2 + g, data = d))), "\n")
pt <- par_template(y ~ x1 + x2 + g, data = d)
cat("frm_simulate identical:",
    identical(q(frm_simulate(y ~ ., d, newparams = pt, seed = 1)),
              q(frm_simulate(y ~ x1 + x2 + g, d, newparams = pt,
                             seed = 1))), "\n")
dnoy <- d[, c("x1", "x2", "g")]
cat("frm_simulate y ~ . with no y column:",
    show(q(dim(frm_simulate(y ~ ., dnoy, newparams = pt, seed = 1)))), "\n")
cat("data = environment:", show(q(frm(y ~ ., data = list2env(d)))), "\n")
cat("data = list:", show(q(class(frm(y ~ ., data = as.list(d))))), "\n")

cat("== update of a dot fit\n")
fit <- frm(y ~ ., data = d)
cat("stored formula:", deparse1(formula(fit)$formula %||% formula(fit)),
    "\n")
d2 <- d; d2$extra <- rnorm(n)
u1 <- q(update(fit, newdata = d2))
cat("update(newdata = d2 with an extra column) fixef:",
    show(if (inherits(u1, "ERR")) u1 else rownames(fixef(u1))), "\n")
u2 <- q(update(fit, . ~ . + I(x1^2)))
cat("update(. ~ . + I(x1^2)) fixef:",
    show(if (inherits(u2, "ERR")) u2 else rownames(fixef(u2))), "\n")
u2b <- q(update(fit, . ~ . - x2))
cat("update(. ~ . - x2) fixef:",
    show(if (inherits(u2b, "ERR")) u2b else rownames(fixef(u2b))), "\n")
dd <- d; dd$yo <- cut(dd$y, c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
dd$y <- NULL
fo <- frm(yo ~ ., data = dd, family = cumulative())
u3 <- q(update(fo, bf(~ ., family = acat())))
cat("update(bf(~ ., family = acat())):",
    show(if (inherits(u3, "ERR")) u3 else
      paste(family(u3)$family, paste(rownames(fixef(u3)), collapse = ","))),
    "\n")
u4 <- q(update(fo, formula. = ~ ., family = acat()))
cat("update(~ ., family = acat()):",
    show(if (inherits(u4, "ERR")) u4 else family(u4)$family), "\n")
u5 <- q(update(fit, bf(~ ., sigma ~ x1)))
cat("update(bf(~ ., sigma ~ x1)):",
    show(if (inherits(u5, "ERR")) u5 else
      paste(rownames(fixef(u5)), collapse = ",")), "\n")
u6 <- q(update(fit, bf(. ~ . + (1 | g))))
cat("update(bf(. ~ . + (1 | g))):",
    show(if (inherits(u6, "ERR")) u6 else
      paste(variables(u6), collapse = ",")), "\n")
# a bf() delta whose family argument and family slot disagree
fg <- frm(y ~ x1, data = d)
u7 <- q(update(fg, bf(~ ., family = student()), family = gaussian()))
cat("bf delta family wins over family =:",
    show(if (inherits(u7, "ERR")) u7 else family(u7)$family), "\n")
# a complete bf() is still a replacement
u8 <- q(update(fg, bf(y ~ x2)))
cat("update(bf(y ~ x2)) fixef:",
    show(if (inherits(u8, "ERR")) u8 else rownames(fixef(u8))), "\n")
# a bf() delta on a model with a dpar formula keeps it
fs <- frm(bf(y ~ x1, sigma ~ x2), data = d)
u9 <- q(update(fs, bf(~ . + x2)))
cat("update(bf(~ . + x2)) of bf(y ~ x1, sigma ~ x2):",
    show(if (inherits(u9, "ERR")) u9 else
      paste(rownames(fixef(u9)), collapse = ",")), "\n")
u10 <- q(update(fs, bf(~ ., sigma ~ 1)))
cat("update(bf(~ ., sigma ~ 1)) replaces sigma:",
    show(if (inherits(u10, "ERR")) u10 else
      paste(rownames(fixef(u10)), collapse = ",")), "\n")
# the delta on a nonlinear model
fnl <- q(frm(bf(y ~ a + b * x1, a ~ 1, b ~ 1, nl = TRUE), data = d))
u11 <- q(update(fnl, bf(~ ., a ~ x2)))
cat("nl update(bf(~ ., a ~ x2)):",
    show(if (inherits(u11, "ERR")) u11 else
      paste(rownames(fixef(u11)), collapse = ",")), "\n")
u12 <- q(update(fnl, bf(~ . + x2)))
cat("nl update(bf(~ . + x2)):", show(if (inherits(u12, "ERR")) u12 else
  deparse1(formula(u12)$formula)), "\n")
# the multivariate refusal
fm <- q(frm(bf(y ~ x1) + bf(x2 ~ x1) + set_rescor(FALSE), data = d))
u13 <- q(update(fm, bf(~ ., sigma ~ x1)))
cat("mv bf delta:", show(if (inherits(u13, "ERR")) u13 else "ok"), "\n")
cat("DONE\n")
