# Lane wt-reunc: how each kind of group-level term stores its
# component, so a re_formula term can be matched against it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(2)
n <- 120
d <- data.frame(g = factor(rep(1:10, length.out = n)),
                h = factor(rep(1:6, each = n / 6)),
                g2 = factor(sample(1:10, n, TRUE)),
                x = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                t = rep(1:12, 10))
d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(6, 0, 0.5)[d$h] + rnorm(n)
A <- diag(10); dimnames(A) <- list(levels(d$g), levels(d$g))
show <- function(lab, f, ...) {
  fit <- tryCatch(suppressWarnings(frm(f, data = d, ...)),
                  error = function(e) {
                    cat(lab, "FIT ERROR", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(invisible())
  cat("==", lab, "\n")
  for (bk in fit$frame$re_blocks) {
    cat(" block", bk$covstruct, "dim", bk$dim, "nlev", length(bk$levels),
        "\n")
    for (cp in bk$components) {
      cat("   comp label=", cp$label, " lp=", cp$lp_key, " off=", cp$offset,
          " dim=", cp$dim, " cnms=", paste(cp$cnms, collapse = ","),
          " bar=", if (is.null(cp$bar)) "NULL" else deparse1(cp$bar),
          " mm=", !is.null(cp$mm),
          if (!is.null(cp$mm)) paste(" mmlabel=", cp$mm$label,
                                     "lhs=", deparse1(cp$mm$lhs)), "\n")
    }
  }
  for (lp in fit$frame$linpreds) {
    for (si in lp$smooths %||% list()) {
      cat("   smooth", si$label, "group_var", si$group_var %||% "-",
          "block_ids", si$block_ids, "\n")
    }
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
show("slope+int", bf(y ~ x + (1 + x | g) + (1 | h)))
show("dbar", bf(y ~ x + (1 + x || g)))
show("inter", bf(y ~ x + (1 | g:h)))
show("nest", bf(y ~ x + (1 | h/g)))
show("id", bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)))
show("gr", bf(y ~ x + (1 | gr(g, cov = A))), data2 = list(A = A))
show("mm", bf(y ~ x + (1 | mm(g, g2))))
show("factor", bf(y ~ x + (1 + f | g)))
show("fs", bf(y ~ s(x, g, bs = "fs", k = 4) + (1 | h)))
show("ar1", bf(y ~ x + ar1(0 + factor(t) | g)))
