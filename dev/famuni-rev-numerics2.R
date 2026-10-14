# Reviewer probe 2: the discrete Weibull mean against an independent
# Euler-Maclaurin sum at M = 1e7 terms; the band of x where an outside
# GEV row's gradient is NaN; GEV mixtures with well-separated clusters.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
dwm <- frmtmb:::dw_moments

cat("\n=== A. dw mean: sum to M = 1e7, then the integral by quadrature\n")
em_ref <- function(mu, k, M = 1e7) {
  L <- log(mu)
  y <- seq_len(M)
  s <- sum(exp(L * y^k))
  c0 <- M + 1
  f <- function(t) exp(L * t^k)
  # the tail integral by quadrature on u = log t, independent of the
  # incomplete-gamma form dw_moments() uses
  ti <- integrate(function(u) exp(u) * f(exp(u)), log(c0), Inf,
                  rel.tol = 1e-13, subdivisions = 2000L)$value
  d0 <- f(c0) * L * k * c0^(k - 1)
  s + ti + f(c0) / 2 - d0 / 12
}
g <- expand.grid(mu = c(0.3, 0.9, 0.99, 0.999), shape = c(0.1, 0.2, 0.35, 0.5))
for (i in seq_len(nrow(g))) {
  m <- dwm(g$mu[i], g$shape[i])$mean
  r <- tryCatch(em_ref(g$mu[i], g$shape[i]), error = function(e) NA)
  cat(sprintf("mu %.3f shape %.2f: dw_moments %.12g  reference %.12g  rel %.2e\n",
              g$mu[i], g$shape[i], m, r, abs(m - r) / abs(r)))
}

cat("\n=== B. outside GEV row, xi 0.06, sigma 1: x where the gradient is NaN\n")
famg <- gen_extreme_value()
xs <- seq(-17, -21, by = -0.01)
p <- c(0, 0, log1p(0.06))
bad <- logical(length(xs)); big <- logical(length(xs))
for (j in seq_along(xs)) {
  y <- xs[j]
  tp <- RTMB::MakeTape(function(p) famg$lpdf(y, list(mu = p[1],
    sigma = exp(p[2]), xi = expm1(p[3]), .eta_sigma = p[2], .eta_xi = p[3]),
    list()), p)
  gr <- tp$jacobian(p)
  bad[j] <- any(!is.finite(gr))
  big[j] <- any(abs(gr) > 1e100)
}
cat("x grid", min(xs), "to", max(xs), "step 0.01:", sum(bad),
    "of", length(xs), "points with a non-finite gradient, at x in",
    paste(range(xs[bad]), collapse = " .. "), ";", sum(big),
    "with |gradient| > 1e100\n")

cat("\n=== C. two-cluster GEV mixtures, clusters 20 sigma apart\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
tab <- NULL
for (s in 1:40) {
  set.seed(7700 + s)
  sep <- 4 + 4 * (s %% 4)
  y <- c(rg(150, 0, 0.3, 0.1), rg(150, sep, 0.3, 0.1))
  w <- character(0)
  f <- tryCatch(withCallingHandlers(
    frm(y ~ 1, family = mixture(gen_extreme_value(), gen_extreme_value()),
        data = data.frame(y = y)),
    warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(f, "error")) {
    tab <- rbind(tab, data.frame(seed = s, sep = sep, conv = NA, ll = NA,
                                 warn = substr(conditionMessage(f), 1, 60)))
  } else {
    tab <- rbind(tab, data.frame(seed = s, sep = sep, conv = f$opt$convergence,
                                 ll = as.numeric(logLik(f)),
                                 warn = if (length(w)) substr(w[1], 1, 60) else ""))
  }
}
print(tab, row.names = FALSE)
cat("errors:", sum(is.na(tab$conv)), " conv != 0:", sum(tab$conv != 0, na.rm = TRUE),
    " warned:", sum(nzchar(tab$warn)), "of", nrow(tab), "\n")
