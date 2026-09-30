# Reviewer probe (lane ceplot): how blocks carry their grouping
# expression (nested g/h, numeric grouping columns), and what
# ce_design_label() returns for them.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
show <- function(fit) {
  for (bk in fit$frame$re_blocks) {
    cat(bk$term_label, "| group_name", bk$group_name, "| bar3",
        deparse(bk$components[[1]]$bar[[3]]), "| levels",
        head(bk$levels, 4), "\n")
  }
}
f1 <- frm(bf(y ~ x + (1 | g / h) + (1 | h)), family = gaussian(), data = dc)
show(f1)
di <- dc
di$g <- as.integer(as.character(di$g))
di$h <- as.integer(as.character(di$h))
f2 <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = di)
show(f2)
ndp <- data.frame(x = 0, g = 1L, h = 1L)
for (k in seq_along(f2$frame$re_blocks)) {
  b <- list(bk = f2$frame$re_blocks[[k]],
            gv = strsplit(f2$frame$re_blocks[[k]]$group_name, ":")[[1]])
  cat("label", k, ":", frmtmb:::ce_design_label(b, ndp), "\n")
}
ndp2 <- data.frame(x = 0, g = 2L, h = 5L)
b <- list(bk = f2$frame$re_blocks[[3]], gv = c("g", "h"))
cat("label g=2,h=5:", frmtmb:::ce_design_label(b, ndp2), "\n")
b <- list(bk = f1$frame$re_blocks[[2]], gv = c("g", "h"))
cat("nested block 2 label at g=1,h=1 factors:",
    frmtmb:::ce_design_label(b, data.frame(x = 0, g = factor(1, 1:6),
                                           h = factor(1, 1:5))), "\n")
# what the design itself calls the row
nd <- data.frame(x = 0, g = 1L, h = 1L)
cat("fitted numeric g=1,h=1 allow_new:\n")
print(tryCatch(fitted(f2, newdata = nd, allow_new_levels = TRUE),
               error = function(e) conditionMessage(e)))
cat("fitted numeric g=1,h=2:\n")
print(fitted(f2, newdata = data.frame(x = 0, g = 1L, h = 2L)))
cat("ce numeric, g=1 h=1 wald:\n")
print(conditional_effects(f2, "x", resolution = 3, re_formula = NULL,
                          conditions = list(g = 1L, h = 1L))$x[, c("g", "h",
                          "estimate__", "lower__", "upper__")])
