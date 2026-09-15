# Punch round 2, nits 1, 4 and 5, source side.
ROOT = "C:/Users/adf44/source/r/frmtmb-wt-generics/"


def edit(rel, pairs, append=None):
    p = ROOT + rel
    s = open(p, encoding="utf-8", newline="").read()
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            raise SystemExit("%s: pattern found %d times: %r"
                             % (rel, n, old[:70]))
        s = s.replace(old, new)
    if append is not None:
        s = s + append
    open(p, "w", encoding="utf-8", newline="\n").write(s)
    print("edited", rel)


# ---- nit 1: the truncation statement is an identity with a formula ---
edit("R/scales.R", [(
"""#' * Under truncation, `fitted()` is the TRUNCATED mean and the
#'   formula is simply wrong. On a `y | trunc(lb = 2000)` lognormal
#'   fit of the same design the naive formula is out by a relative
#'   0.1150 at the worst row (`dev/generics-scale2.R`).
""",
"""#' * Under truncation, `fitted()` is the TRUNCATED mean and the
#'   formula answers a different question. With a lower bound `lb`,
#'   the truncated mean is `exp(mu + sigma^2 / 2) * pnorm(sigma - a) /
#'   pnorm(-a)` with `a = (log(lb) - mu) / sigma`, so the naive
#'   formula falls short by a relative `1 - pnorm(-a) / pnorm(sigma -
#'   a)` in each row. That is an identity, not an estimate, and there
#'   is no single number for it: it depends on where the bound sits in
#'   each row's distribution, near zero for a row far above the bound
#'   and approaching one for a row whose mean is below it
#'   (`dev/generics-trunc.R` sweeps three bounds on one design).
""")])

edit("tests/testthat/test-scale-contract.R", [(
"""test_that("the lognormal identity fails under truncation", {
  # fitted() is the TRUNCATED mean there, so the naive formula is not
  # merely imprecise, it is answering a different question. The
  # assertion is two-sided: the error has to be far above the
  # machine precision the untruncated identity holds to, and the
  # truncated mean has to stay above the truncation bound.
""",
"""test_that("the lognormal identity fails under truncation", {
  # fitted() is the TRUNCATED mean there, so the naive formula is not
  # merely imprecise, it is answering a different question, and how
  # far off it is is itself an identity: per row,
  # (fitted - naive) / fitted = 1 - pnorm(-a) / pnorm(sigma - a) with
  # a = (log(lb) - mu) / sigma. So the assertion is that identity, to
  # machine precision, rather than a chosen threshold on its size,
  # which depends on where the bound sits and has no single value.
"""), (
"""  naive <- exp(mu + sigma(f3)^2 / 2)
  rel <- max(abs(fitted(f3) - naive)) / max(abs(fitted(f3)))
  expect_gt(rel, 1e4 * .Machine$double.eps)
  expect_true(all(fitted(f3) > lb))
""",
"""  s <- sigma(f3)
  fv <- fitted(f3)
  naive <- exp(mu + s^2 / 2)
  a <- (log(lb) - mu) / s
  closed <- 1 - stats::pnorm(-a) / stats::pnorm(s - a)
  expect_lte(max(abs((fv - naive) / fv - closed)),
             8 * .Machine$double.eps)
  # the inverse: the naive formula is not the truncated mean, by far
  # more than the precision the identity above holds to
  expect_gt(max(abs(closed)), 1e6 * .Machine$double.eps)
  expect_true(all(fv > lb))
""")])

# ---- nit 5: the posterior floor ---------------------------------------
edit("DESCRIPTION", [("    posterior,\n", "    posterior (>= 1.0.0),\n")])

# ---- nit 4: the owner-table twin invariant, as a test -----------------
twin = '''
test_that("every borrowed method frmtmb registers has its owner's twin", {
  # A method frmtmb registers in its OWN table on a borrowed name, say
  # `S3method(loo, frmtmb_fit)`, is reachable only while frmtmb's own
  # generic is the one in use. The moment the owner loads, the exported
  # binding is the owner's generic and the method has to be in the
  # OWNER's table too, which is what `S3method(loo::loo, frmtmb_fit)`
  # puts it. A method added with the first directive and not the second
  # passes every test that runs without the owner loaded and goes
  # unreachable the moment it is loaded. Review found the invariant
  # held, 28 rows with 2 deliberate exceptions, and that nothing
  # asserted it (dev/genrev-r2-pairs.R).
  #
  # The two exceptions are frmtmb's own `.default` methods. Registering
  # those on the owner's generic would REPLACE brms's and loo's own
  # defaults for every class, so they stay in frmtmb's table on purpose.
  ns_file <- parseNamespaceFile("frmtmb", dirname(find.package("frmtmb")))
  m <- ns_file$S3methods
  own <- get("frm_generic_owners", envir = asNamespace("frmtmb"))
  missing_twins <- function(m) {
    local <- which(is.na(m[, 4]) & m[, 1] %in% names(own))
    out <- character()
    for (i in local) {
      for (o in own[[m[i, 1]]]) {
        hit <- m[, 1] == m[i, 1] & m[, 2] == m[i, 2] & m[, 4] %in% o
        if (!any(hit)) out <- c(out, paste0(o, "::", m[i, 1], ".", m[i, 2]))
      }
    }
    sort(out)
  }
  deliberate <- c("brms::posterior_summary.default",
                  "loo::loo_compare.default")
  # the guard is only a guard if it has rows to read
  expect_gt(sum(is.na(m[, 4]) & m[, 1] %in% names(own)), 20L)
  expect_equal(missing_twins(m), deliberate)
  # the inverse, built in: drop ONE real twin and the check must name it
  drop <- which(m[, 1] == "loo" & m[, 2] == "frmtmb_fit" & m[, 4] %in% "loo")
  expect_length(drop, 1L)
  expect_equal(missing_twins(m[-drop, , drop = FALSE]),
               sort(c(deliberate, "loo::loo.frmtmb_fit")))
})
'''
edit("tests/testthat/test-generic-collision.R", [], append=twin)
