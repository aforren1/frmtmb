# frm_compat() said "conditional", with a note promising a latent residual, on families
# that have no mu, where frm() refuses (dev/reviews/2026-11-04-release.md, finding 10).

test_that("an autocor term on a family with no mu reads refused, as frm() does", {
  cp <- as.data.frame(frm_compat())
  terms <- c("ar()", "ma()", "arma()", "cosy()", "unstr()")
  fams <- c("rlddm", "bandit2arm_delta", "bandit2arm_dual", "prl_fictitious",
            "bandit4arm2_kalman_filter", "ts_par7", "igt_pvl_delta", "igt_orl")
  for (fm in fams) {
    for (tm in terms) {
      st <- cp$status[(cp$feature_a == tm & cp$feature_b == fm) |
                        (cp$feature_a == fm & cp$feature_b == tm)]
      expect_identical(st, "refused", label = paste(fm, "x", tm))
    }
  }
  # what the row says is what frm() does
  set.seed(7341)
  d <- data.frame(id = rep(1:4, each = 30), trial = rep(1:30, times = 4),
                  choice = stats::rbinom(120, 1, 0.5),
                  pay1 = stats::rbinom(120, 1, 0.7), pay2 = stats::rbinom(120, 1, 0.3))
  expect_error(
    frm(bf(choice | reward(pay1, pay2) ~ ar(time = trial, gr = id, cov = TRUE), tau ~ 1),
        family = bandit2arm_delta(subject = id, trial = trial), data = d),
    "only be specified on 'mu'", fixed = TRUE)
})
