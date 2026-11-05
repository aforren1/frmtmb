# frm_compat() said "conditional", with a note promising a latent residual, on families
# that have no mu, where frm() refuses (dev/reviews/2026-11-04-release.md, finding 10).

test_that("an autocor term on lba and rdm reads refused, as frm() does", {
  cp <- as.data.frame(frm_compat())
  for (fm in c("lba", "rdm")) {
    for (tm in c("ar()", "ma()", "arma()", "cosy()", "unstr()")) {
      st <- cp$status[(cp$feature_a == tm & cp$feature_b == fm) |
                        (cp$feature_a == fm & cp$feature_b == tm)]
      expect_identical(st, "refused", label = paste(fm, "x", tm))
    }
  }
})
