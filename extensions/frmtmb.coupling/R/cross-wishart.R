#' The complex Wishart family for a pair of signals
#'
#' The likelihood of a two-channel cross-spectral matrix. The response
#' is the first auto-spectrum and the rest of the Hermitian matrix rides
#' in `vreal()`, with the degrees of freedom in `vint()`. Its four
#' distributional parameters are the two channel powers, the magnitude
#' squared coherence and the phase, each with its own linear predictor,
#' so a random effect on coherence or a smooth in coherence over
#' frequency is an ordinary [frmtmb::frm()] formula.
#'
#' @return A `frmtmb_family` object.
#'
#' @section The formula this family needs:
#' ```r
#' frm(bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
#'        pow2 ~ 1, coh ~ 1, phase ~ 1),
#'     family = cross_wishart(), data = xs)
#' ```
#' [frm_cross_spectrum()] returns exactly those columns under exactly
#' those names. Every one of the four dpars needs a formula: a dpar left
#' out is a constant, which is right for `pow2` in some designs and
#' almost never right for `coh`.
#'
#' @section The density:
#' For `W` the summed Hermitian matrix at one frequency and `n` its
#' degrees of freedom,
#'
#' ```
#' log f(W) = (n - 2) log det W - tr(S^-1 W) - n log det S - log CGamma_2(n)
#' log CGamma_2(n) = log(pi) + lgamma(n) + lgamma(n - 1)
#' ```
#'
#' In this package's coordinates both hard parts collapse. Writing
#' `S11`, `S22` for the powers, `C` for the coherence and `phi` for the
#' phase, `det S = S11 S22 (1 - C)`, so the log determinant is a sum of
#' three logs and no determinant is taken; and `tr(S^-1 W)` is
#'
#' ```
#' (S22 w11 + S11 w22 - 2 sqrt(C S11 S22) (cos(phi) w12r + sin(phi) w12i))
#'   / (S11 S22 (1 - C))
#' ```
#'
#' which is real arithmetic throughout. **At two channels the complex
#' Wishart needs no complex arithmetic on the tape.** RTMB 1.9's
#' `adcomplex` was used to derive and check this form at general `p`,
#' and it works, including under a Laplace approximation; but at `p = 2`
#' it was measured 5.8 times slower than the closed form for the same
#' answer to six decimals, so the closed form is what ships.
#'
#' @section The trap: a power formula too simple for the data:
#' The likelihood couples the powers and the coherence. A power model too
#' simple for the spectrum does not merely leave `mu` and `pow2` wrong:
#' it reweights the rows each coherence parameter pools, and moves the
#' coherence with them. It converges without complaint.
#'
#' **Check it by comparing AIC with and without a frequency term on the
#' powers.** Nothing here does it for you. A draft of this package
#' warned automatically from a `fit_check` hook that measured the spread
#' of the power residuals, and that statistic is not monotone in the
#' coherence displacement it was meant to catch: it went silent at a
#' 2.7-fold power spread where the reported coherence was already nearly
#' double the truth, and it fired at a ratio of 20.5 on a correctly
#' specified model whose coherences were right. It was withdrawn rather
#' than retuned, because a check that punishes correct models teaches
#' users to ignore it.
#'
#' An AIC comparison does work, and works where the hook did not. Two
#' bands over 400 frequencies at `n = 8`, true coherences 0.35 and 0.05,
#' `coh ~ band`, varying the power spread across the record:
#'
#' | power spread | flat powers report | correct powers report | delta AIC |
#' | --- | --- | --- | --- |
#' | 1.6x | 0.295, 0.083 | 0.330, 0.068 | 80 |
#' | 2.7x | 0.268, 0.078 | 0.342, 0.041 | 487 |
#' | 4.4x | 0.256, 0.136 | 0.344, 0.054 | 954 |
#' | 7.3x | 0.240, 0.206 | 0.342, 0.047 | 1849 |
#' | 19.6x | 0.229, 0.458 | 0.341, 0.052 | 3716 |
#' | 385x | 0.222, 0.875 | 0.342, 0.059 | 11518 |
#'
#' The flat fit reports the high band as nine times its true coherence
#' at a 19.6-fold spread and as seventeen times at 385-fold, and AIC
#' rejects it by 80 or more even at the mildest spread tested. The
#' correctly specified fits recover 0.35 and 0.05 throughout.
#'
#' So: fit the powers twice, once with a frequency term and once
#' without, and read the coherence off the one AIC prefers.
#' **Model the powers before reading any coherence.**
#'
#' @section Why the links are the constraint, in floating point too:
#' The spectral matrix must be positive definite at every frequency, in
#' every group, for every draw the optimizer tries. `det S > 0` needs
#' only `C < 1`, and `C` reaches the density through a logit link, so no
#' value of any linear predictor can produce an invalid matrix **in
#' exact arithmetic**. No eigenvalue floor and no projection step are
#' applied, and none is needed.
#'
#' Exact arithmetic is not what runs. `plogis(eta)` is exactly 1 in
#' double precision from about `eta = 36.74` upward, so a density that
#' formed `1 - C` by subtraction would read `log(0)` and divide by zero
#' there, and would already be wrong in the third digit at `eta = 30`.
#' That region is reachable: a random effect on `coh` with one group
#' near a true coherence of 1 drives the estimate past `eta = 34`.
#'
#' This family therefore never subtracts. Core stores each dpar's linear
#' predictor beside it, and the logit link's `logit_eta` field makes the
#' log-odds exact, so `log(1 - C)` is computed as
#' `-logspace_add(0, eta)` and `1 / (1 - C)` as its exponential. Both
#' are exact to `eta = 709`, where the double range itself ends, against
#' `36.74` for the naive form. Measured against a high-precision
#' reference at `n = 16`: relative error 1.0e-03 at `eta = 30` and NaN
#' at `eta = 40` before the change, below 1e-15 at both after it.
#'
#' One path keeps the plain round trip, and says so rather than being
#' floored: `residuals(type = "deviance")` runs off the tape, where core
#' does not store the linear predictor, so it uses `log1p(-C)`. That is
#' accurate until `C` rounds to exactly 1, which needs `eta` past 36.74
#' and is further than a fit reaches in practice. Measured on the
#' degenerate case above, two signals differing by 1e-5 of noise: the
#' fit lands at `eta = 23.0`, `C` is 0.99999999989743915 rather than 1,
#' and all 127 response, Pearson and deviance residuals are finite.
#' Past 36.74 they would be `NaN`, and no floor is applied to hide it.
#'
#' @section What it refuses:
#' * More than two channels. This is a scope decision, and **core is not
#'   the obstacle**: `frm()` carries a matrix-valued response, and both
#'   a gaussian fit and a custom family indexing `y[, 1]` and `y[, 2]`
#'   work today. What is missing is here rather than there. A
#'   `p`-channel model needs `p^2` linear predictors written out by
#'   hand; the clean `power, power, coherence, phase` coordinates of
#'   this family do not survive past `p = 2`, so it would be written on
#'   raw Cholesky entries that mean nothing to a reader; and the log
#'   determinant and the trace stop being closed forms, which puts
#'   `adcomplex` back on the tape at about six times the cost. The
#'   density itself generalizes and was checked at `p = 4`. Fit pairs.
#' * `n < 2`. A single complex draw gives a rank-one matrix whose
#'   coherence is exactly 1 and whose log determinant is `-Inf`; the
#'   normalizing constant is `Inf` there as well, since
#'   `lgamma(n - 1)` reaches a non-positive argument. This is checked
#'   in `valid_y()`, before any tape is built.
#' * `simulate()`. A draw from this family is a whole Hermitian matrix
#'   and the response slot carries one number of the four, so a
#'   simulated `w11` with the observed `w22`, `w12r` and `w12i` beside
#'   it is not a draw from anything. [frm_cross_simulate()] returns all
#'   four columns instead.
#'
#' @seealso [frm_cross_spectrum()] to build the response,
#'   [frm_coherence()] and [frm_phase()] to read the fit
#' @examples
#' set.seed(4)
#' src <- rnorm(4096)
#' xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
#'                          segments = 16)
#' xs <- xs[xs$freq < 0.1, ]
#' fit <- frmtmb::frm(
#'   frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
#'              pow2 ~ 1, coh ~ 1, phase ~ 1),
#'   family = cross_wishart(), data = xs)
#' frm_coherence(fit)[1, ]
#' @export
cross_wishart <- function() {
  custom_family(
    "cross_wishart",
    dpars = c("mu", "pow2", "coh", "phase"),
    links = list(mu = "log", pow2 = "log", coh = "logit",
                 phase = "identity"),
    lpdf = cw_lpdf,
    valid_y = cw_valid_y,
    init_dpars = list(mu = function(y, aterms) mean(y / aterms$vint1),
                      pow2 = function(y, aterms)
                        mean(aterms$vreal1 / aterms$vint1),
                      coh = cw_init_coh,
                      phase = function(y, aterms)
                        atan2(sum(aterms$vreal3), sum(aterms$vreal2))),
    type = "continuous",
    post = list(mean_fn = function(dpars, aterms) aterms$vint1 * dpars$mu,
                var_fn = function(dpars, aterms) aterms$vint1 * dpars$mu^2,
                dev_fn = cw_dev),
    sim_refusal = paste0(
      "a draw is a whole Hermitian matrix and the response carries only ",
      "its first entry; frmtmb.coupling::frm_cross_simulate() returns ",
      "all four columns"),
    required_aterms = c("vint1", "vreal1", "vreal2", "vreal3"),
    accepts_aterms = c("vint", "vreal", "weights"))
}

#' The coherence complement, `1 - C`, without ever subtracting from 1.
#'
#' `C` reaches a density through the logit link's `linkinv`, and
#' `plogis(eta)` is exactly 1 in double precision from about
#' `eta = 36.74` upward. Forming `1 - C` there gives exactly 0, so
#' `log(1 - C)` is `-Inf` and the trace divides by zero: the density and
#' its gradient are both NaN where the truth is an ordinary large
#' negative number. Below that the subtraction is merely inaccurate,
#' 1e-3 relative at `eta = 30`.
#'
#' Core anticipated this. `build_objective()` stores each dpar's linear
#' predictor beside it as `.eta_<dpar>`, and the logit link's
#' `logit_eta` field is the identity, so the log-odds are available
#' exactly on the tape. `log(1 - plogis(eta))` is then
#' `-logspace_add(0, eta)`, which is exact everywhere, and `1 / (1 - C)`
#' is its exponential, finite to `eta = 709` rather than to 36.74.
#'
#' Off the tape the `.eta_` entries are absent by design and the plain
#' round trip is used. `cw_dev()` is the only caller that runs there,
#' and `log1p(-C)` holds until `C` rounds to exactly 1 above
#' `eta = 36.74`, which is past where a fit lands: a deliberately
#' degenerate fit reached `eta = 23.0` with every residual finite.
#'
#' @noRd
cw_complement <- function(dpars) {
  eta <- dpars[[".eta_coh"]]
  if (is.null(eta)) {
    ch <- dpars$coh
    return(list(log = log1p(-ch), inv = 1 / (1 - ch)))
  }
  lg <- -RTMB::logspace_add(0 * eta, eta)
  list(log = lg, inv = exp(-lg))
}

#' The log density, vectorized over rows and written for the tape.
#'
#' Everything here is real arithmetic; see the family's "The density"
#' section for why that is enough at two channels. The one place that
#' needs care is `1 - C`, which is never formed by subtraction; see
#' `cw_complement()` above.
#'
#' @noRd
cw_lpdf <- function(y, dpars, aterms) {
  n <- aterms$vint1
  w22 <- aterms$vreal1; w12r <- aterms$vreal2; w12i <- aterms$vreal3
  s11 <- dpars$mu; s22 <- dpars$pow2; ch <- dpars$coh; ph <- dpars$phase
  cmp <- cw_complement(dpars)
  logdetW <- log(y * w22 - w12r^2 - w12i^2)
  logdetS <- log(s11) + log(s22) + cmp$log
  recross <- sqrt(ch * s11 * s22) * (cos(ph) * w12r + sin(ph) * w12i)
  tr <- (s22 * y + s11 * w22 - 2 * recross) * cmp$inv / (s11 * s22)
  (n - 2) * logdetW - tr - n * logdetS -
    (log(pi) + lgamma(n) + lgamma(n - 1))
}

#' The unit deviance, against the saturated fit S = W / n.
#'
#' At the saturated value `tr(S^-1 W) = 2n` and
#' `log det S = log det W - 2 log n`, so the whole expression is
#' available in closed form and a deviance residual means what it means
#' for any other family.
#'
#' @noRd
cw_dev <- function(y, dpars, aterms) {
  n <- aterms$vint1
  w22 <- aterms$vreal1; w12r <- aterms$vreal2; w12i <- aterms$vreal3
  s11 <- dpars$mu; s22 <- dpars$pow2; ch <- dpars$coh; ph <- dpars$phase
  cmp <- cw_complement(dpars)
  logdetW <- log(y * w22 - w12r^2 - w12i^2)
  logdetS <- log(s11) + log(s22) + cmp$log
  recross <- sqrt(ch * s11 * s22) * (cos(ph) * w12r + sin(ph) * w12i)
  tr <- (s22 * y + s11 * w22 - 2 * recross) * cmp$inv / (s11 * s22)
  2 * (tr - 2 * n + n * (logdetS - logdetW + 2 * log(n)))
}

#' The moment start for coherence, pooled over rows.
#'
#' A per-row coherence is exactly 1 whenever n is 2, so the start is
#' taken from the summed matrix and then pulled off both boundaries.
#'
#' @noRd
cw_init_coh <- function(y, aterms) {
  num <- sum(aterms$vreal2)^2 + sum(aterms$vreal3)^2
  den <- sum(y) * sum(aterms$vreal1)
  min(max(num / den, 0.02), 0.9)
}

#' @noRd
cw_valid_y <- function(y, aterms) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("the response of cross_wishart() is an auto-spectrum and must ",
         "be positive and finite in every row.", call. = FALSE)
  }
  n <- aterms$vint1
  w22 <- aterms$vreal1; w12r <- aterms$vreal2; w12i <- aterms$vreal3
  if (any(!is.finite(w22)) || any(w22 <= 0)) {
    stop("the second auto-spectrum, vreal1, must be positive and finite ",
         "in every row.", call. = FALSE)
  }
  if (any(!is.finite(w12r)) || any(!is.finite(w12i))) {
    stop("the cross-spectrum, vreal2 and vreal3, holds a missing or ",
         "non-finite value.", call. = FALSE)
  }
  bad <- which(n < 2)
  if (length(bad)) {
    stop("cross_wishart() needs at least 2 degrees of freedom per row ",
         "and row ", bad[1L], " has ", n[bad[1L]],
         ". One complex draw gives a rank-one matrix whose coherence is ",
         "exactly 1 whatever the signals did; use more segments, more ",
         "tapers, or a wider smooth.", call. = FALSE)
  }
  det <- y * w22 - w12r^2 - w12i^2
  bad <- which(!(det > 0))
  if (length(bad)) {
    stop("the cross-spectral matrix in row ", bad[1L],
         " is not positive definite: its determinant is ",
         signif(det[bad[1L]], 4),
         ". A coherence of 1 to machine precision does this, and so does ",
         "assembling the four columns in the wrong order.", call. = FALSE)
  }
  if (any(n < 4)) {
    warning("cross_wishart() has fewer than 4 degrees of freedom in ",
            sum(n < 4), " of ", length(n),
            " rows. The density exists at 2, but a coherence from that ",
            "few draws carries a standard deviation near 0.29 and its ",
            "variance component will not identify.", call. = FALSE)
  }
  invisible(TRUE)
}
