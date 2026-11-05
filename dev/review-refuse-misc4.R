# Item 2, 3, 7: (d) and (e) of the decision; brms on positive mixtures and on
# mi() x cens() by brm(); the thres mv reason. Seed 10500.
source("dev/review-refuse-lib.R")

# (d) a response read on some rows only. z itself missing on 40 of the 60
# rows where y is missing, so z's density never reads those rows' y
set.seed(9503)
n <- 200
x <- rnorm(n); eta <- 0.3 + 0.5 * x
d <- data.frame(x = x, y = eta + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
rows <- sort(sample(n, 60))
d$y[rows] <- NA
d$z[rows[1:40]] <- NA
j <- tr(frm(bf(y | mi() ~ x) + exgaussian() + bf(z | mi() ~ mi(y)) + gaussian(), data = d))
d2 <- d[-rows[1:40], ]
j2 <- tr(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) + gaussian(), data = d2))
pr("(d) %s | %s", msg(j, 150), msg(j2, 150))
if (!is_err(j) && !is_err(j2)) {
  ej <- unlist(j$estimates[c("beta", "betad")]); e2 <- unlist(j2$estimates[c("beta", "betad")])
  pr("(d) read on 20 of 60 missing rows; logLik %.3f vs %.3f, gap per dropped row %+.4f", as.numeric(logLik(j)), as.numeric(logLik(j2)),
     (as.numeric(logLik(j)) - as.numeric(logLik(j2))) / 40)
  print(round(rbind(all_rows = ej, unread_dropped = e2), 3))
}
# (e) student at 30 percent
set.seed(9504)
x <- rnorm(200); d <- data.frame(x = x, y = 0.3 + 0.5 * x + rt(200, 3)); rows <- sample(200, 60); d$y[rows] <- NA
a <- frm(bf(y | mi() ~ x), family = student(), data = d)
b <- frm(bf(y ~ x), family = student(), data = d[-rows, ])
print(names(a$estimates$betad))
print(unlist(a$estimates[c("beta", "betad")])); print(unlist(b$estimates[c("beta", "betad")]))
vc <- tr(as.data.frame(frmtmb::fixef(a)))
pr("(e) nu as summary reads it: joint %s; observed rows %s", format(tr(frmtmb:::fit_dpars <- NULL), digits = 3), "")

# brms: a positive mixture with mi(): bound on Ymi
set.seed(10500)
n <- 80; x <- rnorm(n); kk <- rbinom(n, 1, 0.4)
d2 <- data.frame(x = x, y = ifelse(kk == 1, rlnorm(n, 1, 0.4), rlnorm(n, -1, 0.4)))
d2$y[seq(5, n, 15)] <- NA
for (fx in list(brms::mixture(brms::lognormal(), brms::lognormal()), brms::mixture(brms::Beta(), brms::Beta()))) {
  dd <- d2
  if (identical(fx$mix[[1]]$family, "beta")) dd$y <- plogis(dd$y - 1)
  b <- tr(brms::make_stancode(brms::bf(y | mi() ~ x), family = fx, data = dd))
  pr("brms mixture mi: %s; %s", msg(b, 100), if (is_err(b)) "" else paste(trimws(grep("Ymi;", strsplit(b, "\n")[[1]], value = TRUE)), collapse = " "))
}
# brms bounds for gev, exgaussian, skew_normal, asym_laplace, student
for (ft in c("gen_extreme_value", "exgaussian", "skew_normal", "asym_laplace", "student", "von_mises", "beta", "frechet")) {
  dd <- d2; if (ft == "beta") dd$y <- plogis(dd$y - 1); if (ft == "von_mises") dd$y <- atan(dd$y)
  b <- tr(brms::make_stancode(brms::bf(y | mi() ~ x), family = brms::brmsfamily(ft), data = dd))
  pr("brms %-20s %s", ft, if (is_err(b)) msg(b, 100) else paste(trimws(grep("Ymi;", strsplit(b, "\n")[[1]], value = TRUE)), collapse = " "))
}

# brm() on mi() x cens()
set.seed(10004); n <- 180; x <- rnorm(n)
d10 <- data.frame(x = x, y = rnorm(n, 0.5 * x), ce = 0L); d10$y[seq(5, n, 9)] <- NA
b1 <- tr(suppressWarnings(brms::brm(brms::bf(y | mi() + cens(ce) ~ x), data = d10, chains = 1, iter = 100, refresh = 0, seed = 1, silent = 2, init = 0)))
pr("brm() mi() + cens(): %s; draws %s", msg(b1, 200), if (is_err(b1)) "-" else tryCatch(brms::ndraws(b1), error = function(e) conditionMessage(e)))
# the same with the missing rows marked right- or left-censored
for (cc in c(1L, -1L)) {
  dd <- d10; dd$ce[is.na(dd$y)] <- cc
  sd0 <- brms::make_standata(brms::bf(y | mi() + cens(ce) ~ x), data = dd)
  pr("mi + cens, missing rows coded %+d: Jevent has them %s; Jrcens %s; Jlcens %s; Y there %s", cc,
     any(which(is.na(dd$y)) %in% sd0$Jevent), any(which(is.na(dd$y)) %in% sd0$Jrcens), any(which(is.na(dd$y)) %in% sd0$Jlcens), unique(sd0$Y[is.na(dd$y)]))
}

# thres x mvbf with a shared key: read the survey's reason
set.seed(10003); n <- 180; x <- rnorm(n); g <- factor(sample(1:10, n, TRUE))
o <- cut(0.7 * x + rnorm(n), 4, labels = FALSE)
d9 <- data.frame(x = x, o = o, w = rnorm(n, 1 + 0.3 * x), g = g)
E$eq_specs[["variant"]] <- function() list(bf = "mvbf(bf(o | thres(3) ~ x + (1 | p | g), family = cumulative()), bf(w ~ x + (1 | p | g)), rescor = FALSE)", family = NA, data = d9)
r <- E$eq_measure("variant", list(status = "untested"))
pr("thres + key mv: verdict %s; reason %s; l4 %s", r$verdict, substr(r$reason %||% "", 1, 200), paste(format(unlist(r$l4)[seq_len(min(6, length(unlist(r$l4))))], digits = 3), collapse = ","))
