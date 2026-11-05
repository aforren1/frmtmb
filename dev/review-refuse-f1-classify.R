# Short check, item 1 and 3: try to fool mi_rows_read(). Each design says how
# many missing values a density really reads; the fit says how many it held.
# Then the placeholder is moved far inside the support: a logLik that moves
# means it reached the likelihood. Reference where one exists: the same model
# with the unread rows deleted by hand. Seed 11000.
source("dev/review-refuse-lib.R")
suppressPackageStartupMessages(library(frmtmb.spatial))
set.seed(11000)
n <- 120
x <- rnorm(n); g <- factor(rep(1:12, 10)); t0 <- rep(1:10, 12)
y0 <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))
d <- data.frame(x = x, g = g, t0 = t0, f = factor(sample(1:3, n, TRUE)), ord = sample(1:4, n, TRUE),
                y = y0, yn = 0.3 + 0.5 * x + rnorm(n), s = 0.3, id = seq_len(n))
d$z <- 0.5 + 0.6 * y0 + 0.3 * x + rnorm(n, 0, 0.7)
d$w <- rnorm(n, 0.4 * x)
M <- seq(4, n, by = 12)        # 10 rows missing in y
d$y[M] <- NA; d$yn[M] <- NA
fit_one <- function(lab, expr, want, ref = NULL) {
  ori <- get("mi_placeholder", asNamespace("frmtmb"))
  a <- tr(expr())
  if (is_err(a)) { pr("%-34s refused: %s", lab, msg(a, 120)); return(invisible(NULL)) }
  nl <- length(a$estimates$miss %||% numeric())
  assignInNamespace("mi_placeholder", function(obs) stats::quantile(obs, 0.93, names = FALSE), "frmtmb")
  b <- tr(expr())
  assignInNamespace("mi_placeholder", ori, "frmtmb")
  dll <- if (is_err(b)) NA else as.numeric(logLik(a)) - as.numeric(logLik(b))
  ec <- if (!is.null(ref)) {
    r <- tr(ref())
    if (is_err(r)) msg(r, 80) else sprintf("%.3e", as.numeric(logLik(a)) - as.numeric(logLik(r)))
  } else "-"
  pr("%-34s want %2s latent, holds %2d; nobs %d/%s; logLik %.4f; placeholder moves it by %.1e; minus hand deletion %s", lab, format(want), nl,
     nobs(a), format(attr(logLik(a), "nobs")), as.numeric(logLik(a)), dll, ec)
  invisible(a)
}
Gam <- stats::Gamma(link = "log")
DM <- d[-M, ]

# the plain cases
fit_one("unread, Gamma", function() frm(bf(y | mi() ~ x), family = Gam, data = d), 0, function() frm(bf(y ~ x), family = Gam, data = DM))
fit_one("read by z", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y) + x) + gaussian(), data = d), 10)
fit_one("read by dpar sigma of z", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ x, sigma ~ mi(y)) + gaussian(), data = d), 10)
fit_one("read in mi(y):x", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y):x + x) + gaussian(), data = d), 10)
fit_one("read in mi(y) * x", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y) * x) + gaussian(), data = d), 10)
fit_one("read in nlpar", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ a + b, a ~ mi(y), b ~ 0 + x, nl = TRUE) + gaussian(), data = d), 10)
fit_one("read in mixture mu1", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ 1, mu1 ~ mi(y), mu2 ~ 1) + mixture(gaussian(), gaussian(), order = "none"), data = d), 10)
fit_one("read in random slope", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ x + (0 + mi(y) | g)) + gaussian(), data = d), 10)
fit_one("read by mi(y) in s()", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ s(mi(y))) + gaussian(), data = d), 10)
fit_one("read by I(mi(y)^2)", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ I(mi(y)^2)) + gaussian(), data = d), 10)
fit_one("read by mo(ord):mi(y)", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mo(ord):mi(y)) + gaussian(), data = transform(d, ord = ordered(ord))), 10)
fit_one("read by cs(mi(y))", function() frm(bf(y | mi() ~ x) + Gam + bf(ord ~ cs(mi(y))) + cumulative(), data = d), 10)
fit_one("read twice", function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y)) + gaussian() + bf(w ~ mi(y)) + gaussian(), data = d), 10)
fit_one("mi(y) by itself in its dpar", function() frm(bf(y | mi() ~ x, shape ~ mi(y)) + Gam, data = d), 10)
fit_one("nl body names y", function() frm(bf(z ~ b0 + b1 * y, b0 ~ 1, b1 ~ 1, nl = TRUE) + gaussian() + bf(y | mi() ~ x) + Gam, data = d), "?")
fit_one("mi(s), some NA", function() frm(bf(y | mi(s) ~ x), family = Gam, data = d), n)

# plain uses of the column, which must not see a placeholder
for (lab in c("z ~ y plain", "z | se(y) ~ x", "z | weights(y) ~ x", "z | subset(ok) ~ x")) {
  fm <- switch(lab,
    `z ~ y plain` = bf(z ~ y + x),
    `z | se(y) ~ x` = bf(z | se(y) ~ x),
    `z | weights(y) ~ x` = bf(z | weights(y) ~ x),
    `z | subset(ok) ~ x` = bf(z | subset(ok) ~ x))
  dd <- transform(d, ok = !is.na(y))
  r0 <- tr(frm(fm, family = gaussian(), data = dd[!is.na(dd$y), ]))
  fj <- function() frm(bf(y | mi() ~ x) + Gam + fm + gaussian(), data = dd)
  a <- fit_one(paste("beside mv:", lab), fj, 0)
  if (!is.null(a) && !is_err(r0)) {
    y1 <- frm(bf(y ~ x), family = Gam, data = dd[!is.na(dd$y), ])
    pr("   sum of the two separate complete-case fits %.4f (z on %d rows)", as.numeric(logLik(r0)) + as.numeric(logLik(y1)), nobs(r0))
  }
}
ba <- tr(brms::make_standata(brms::bf(y | mi() ~ x, family = Gam) + brms::bf(z ~ y + x) + brms::set_rescor(FALSE), data = d))
pr("brms, z ~ y plain beside mi(y): %s; N_z %s", msg(ba, 110), if (is_err(ba)) "-" else ba$N_z)

# a chain: z itself missing where y is. z's density is taken out, so y's row is
# read by nothing
dz <- d; dz$z[M[1:6]] <- NA
fit_one("chain: z missing on 6 of 10", function() frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y) + x) + gaussian(), data = dz), 4,
        function() frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y) + x) + gaussian(), data = dz[-M[1:6], ]))
# subset and idx: mi(y, idx = id2) reads rows through a map
dd <- d; dd$ok <- TRUE; dd$ok[M[1:6]] <- FALSE; dd$z2 <- ifelse(dd$ok, dd$z, 0)
sub <- dd[dd$ok, ]
fit_one("idx: read on 4 rows of 10", function() frm(bf(y | mi() + index(id) ~ x) + Gam + bf(z2 | subset(ok) ~ mi(y, idx = id) + x) + gaussian(), data = dd), 4,
        function() {
          d3 <- dd[dd$ok, ]
          frm(bf(y | mi() ~ x) + Gam + bf(z2 ~ mi(y) + x) + gaussian(), data = d3)
        })
# rows a weight of 0 turns off in the reader; weights on the mi response
fit_one("reader weight 0 on missing rows", function() frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ mi(y) + x) + gaussian(), data = transform(d, wz = ifelse(is.na(y), 0, 1))), "0?")
fit_one("weights on the unread", function() frm(bf(y | mi() + weights(wy) ~ x), family = Gam, data = transform(d, wy = 1 + (x > 0))), 0,
        function() frm(bf(y | weights(wy) ~ x), family = Gam, data = transform(d, wy = 1 + (x > 0))[-M, ]))

# autocorrelation: the row's value enters a joint density or a lagged residual
od <- order(d$g, d$t0); da <- d[od, ]
lab_ac <- list(
  "ar(cov = FALSE)" = "bf(yn | mi() ~ x + ar(time = t0, gr = g, p = 1))",
  "ar(cov = TRUE)" = "bf(yn | mi() ~ x + ar(time = t0, gr = g, p = 1, cov = TRUE))",
  "ma(cov = TRUE)" = "bf(yn | mi() ~ x + ma(time = t0, gr = g, q = 1, cov = TRUE))",
  "arma(cov = TRUE)" = "bf(yn | mi() ~ x + arma(time = t0, gr = g, cov = TRUE))",
  "cosy" = "bf(yn | mi() ~ x + cosy(time = t0, gr = g))",
  "unstr" = "bf(yn | mi() ~ x + unstr(time = t0, gr = g))")
for (nm in names(lab_ac)) {
  fm <- eval(parse(text = lab_ac[[nm]]), envir = asNamespace("frmtmb"))
  a <- fit_one(paste("autocor", nm), function() frm(fm, data = da), 10)
  b <- tr(brms::make_stancode(eval(parse(text = lab_ac[[nm]]), envir = asNamespace("brms")), data = da))
  pr("   brms builds: %s; %s", msg(b, 90), if (is_err(b)) "" else paste(trimws(grep("Yl|Ymi", strsplit(b, "\n")[[1]], value = TRUE))[1:2], collapse = " | "))
}
W <- 0.7^abs(outer(1:n, 1:n, "-")); M2 <- W
a <- fit_one("fcor", function() frm(bf(yn | mi() ~ x + fcor(M2)), data = da, data2 = list(M2 = M2)), 10)
A <- (abs(outer(seq_len(n), seq_len(n), "-")) == 1) * 1
rownames(A) <- colnames(A) <- NULL
a <- fit_one("sar", function() frm(bf(yn | mi() ~ x + sar(A)), data = da, data2 = list(A = A)), 10)

# rescor: refused or latent
fit_one("rescor, y latent", function() frm(mvbf(bf(yn | mi() ~ x), bf(w ~ x), rescor = TRUE), data = d), 10)
fit_one("rescor among readers", function() frm(bf(yn | mi() ~ x) + gaussian() + mvbf(bf(z ~ mi(yn)), bf(w ~ x), rescor = TRUE), data = d), "?")
# an hmm or group mixture holds the row in a joint density
fit_one("mixture(groups = ) with mi", function() frm(bf(yn | mi() ~ x), family = mixture(gaussian(), gaussian(), groups = ~ g), data = d), 10)
