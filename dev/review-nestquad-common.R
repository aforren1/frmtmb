# Shared setup for the nestquad review scripts. Source with ARM set to
# "lane" (default) or "base" (rellib-r18 only).
arm <- Sys.getenv("ARM", "lane")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad"
lane_lib <- Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/wt-nestquad-lib")
base_libs <- c("C:/Users/adf44/source/r/rellib-r18",
               "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "base") base_libs else c(lane_lib, base_libs))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", arm, format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
has_nest <- exists("quad_nest_fit", ns)
cat("has quad_nest_fit:", has_nest, "\n")

# peak working set of this process, MB (ps), for memory checks
peak_mb <- function() {
  m <- ps::ps_memory_info()
  unname(m[["peak_wset"]] %||% m[["rss"]]) / 2^20
}
`%||%` <- function(a, b) if (is.null(a)) b else a
cur_mb <- function() unname(ps::ps_memory_info()[["rss"]]) / 2^20

lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

# Gauss-Hermite for weight exp(-x^2) via statmod-free Golub-Welsch,
# physicists' convention, written independently of the package
gh_phys <- function(k) {
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

# log int exp(f(v)) dv by adaptive GH: center at the mode of f (found by
# optimize on a wide bracket), scale by the numerical curvature
adapt_int <- function(f, k, lo = -30, hi = 30) {
  m <- optimize(f, c(lo, hi), maximum = TRUE, tol = 1e-10)$maximum
  h <- 1e-4
  c2 <- -(f(m + h) - 2 * f(m) + f(m - h)) / h^2
  if (!is.finite(c2) || c2 <= 0) c2 <- 1
  sc <- sqrt(2 / c2)
  g <- gh_phys(k)
  v <- m + sc * g$x
  fv <- vapply(v, f, 1)
  log(sc) + lse(log(g$w) + g$x^2 + fv)
}

# Nested reference from the model alone. eta0: fixed linear predictor
# per row; za, zi: design values of the outer and inner effect per row;
# ga, gi: outer and inner group of each row (inner labels unique);
# lpdf(eta, rows): weighted log density of `rows` at predictor `eta`.
nest_ref <- function(eta0, za, zi, ga, gi, sd_a, sd_i, lpdf, k = 40) {
  tot <- 0
  for (a in unique(ga)) {
    ra <- which(ga == a)
    inn <- unique(gi[ra])
    fa <- function(u) {
      s <- dnorm(u, 0, sd_a, log = TRUE)
      for (c in inn) {
        rc <- ra[gi[ra] == c]
        fc <- function(v) {
          dnorm(v, 0, sd_i, log = TRUE) +
            lpdf(eta0[rc] + za[rc] * u + zi[rc] * v, rc)
        }
        s <- s + adapt_int(fc, k, -10 * sd_i - 10, 10 * sd_i + 10)
      }
      s
    }
    tot <- tot + adapt_int(fa, k, -10 * sd_a - 10, 10 * sd_a + 10)
  }
  tot
}

# block sds at outer parameters par (the fit's own covstruct map)
sds_of <- function(fit, par = fit$opt$par) {
  th <- ns$imp_par_list(fit$frame$par_template, par, fit$frame$map)$theta
  s <- vapply(fit$frame$re_blocks, function(bk) {
    sqrt(ns$cs_entry(bk$covstruct)$vcov(th[bk$theta_idx], bk)[1, 1])
  }, 1)
  names(s) <- vapply(fit$frame$re_blocks, `[[`, "", "term_label")
  s
}

# The lane's nested objective at par with given nodes, anchored at par
nest_at <- function(fit, par, nodes) {
  fr <- fit$frame
  lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template,
                         random = "b", map = fr$map, silent = TRUE)
  plan <- ns$quad_check_blocks(fr)
  anc <- ns$quad_nest_anchor(lap, plan, par)
  pl <- ns$imp_par_list(fr$par_template, par, fr$map)
  qo <- ns$quad_nest_objective(fr, plan, anc, nodes, pl)
  qo$prime(pl)
  -qo$fn(pl)
}
lap_at <- function(fit, par) {
  fr <- fit$frame
  lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template,
                         random = "b", map = fr$map, silent = TRUE)
  -lap$fn(par)
}
try_msg <- function(expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(x) {
      w <<- c(w, paste(class(x)[1], conditionMessage(x)))
      invokeRestart("muffleWarning")
    })
  attr(r, "warns") <- w
  r
}
