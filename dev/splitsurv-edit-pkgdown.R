# Core's _pkgdown.yml for the moves: the curve pages and article into
# core, cox_baseline out, the survival subsite in the menu, and
# redirects from the pages that moved. Run once.
f <- "_pkgdown.yml"
x <- readLines(f)
one <- function(s) {
  i <- which(x == s)
  if (length(i) != 1L) stop("not once: ", s)
  i
}
# the menu: frmtmb.survival, after frmtmb.sample
i <- one("      - text: frmtmb.sample (NUTS sampling and posteriors)")
x <- append(x, c("      - text: frmtmb.survival (survival families)",
                 "        href: https://aforren1.github.io/frmtmb/frmtmb.survival/"),
            after = i + 1L)
# cox_baseline() is frmtmb.survival's now
x <- x[-one("  - cox_baseline")]
# the curve pages, beside the seam they read
i <- one("  - frm_extra_cov_deriv")
x <- append(x, c("  - frm_curve", "  - frm_curve_deriv", "  - frm_curve_feature"),
            after = i)
# the article
i <- one("  - case-studies")
x <- append(x, "  - curve-inference", after = i)
# redirects: the core page that moved out, and the moved pages of the
# subsite the functions came from
i <- one("redirects:")
x <- append(x, c(
  "  - [\"reference/cox_baseline.html\", \"frmtmb.survival/reference/cox_baseline.html\"]"),
  after = i)
writeLines(x, f)
