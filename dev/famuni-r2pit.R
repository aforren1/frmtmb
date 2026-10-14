p <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/rtmb-pitfalls.md"
x <- paste(readLines(p), collapse = "\n")
old <- "    with no constant operand: with `u = x + abs(x)`, `u / u - 1` is 0
    with derivatives 0 for `x > 0`, NaN for `x <= 0`, and silent on
    doubles (core's `nan_unless_pos()`)."
new <- "    with no constant operand: with `u = x + abs(x)`, `u / u - 1` is 0
    for `x > 0`, NaN for `x <= 0`, and silent on doubles (core's
    `nan_unless_pos()`). Its derivatives are 0 where it is 0 except at
    the ends of the doubles: the gradient is NaN for `0 < x < 2.3e-308`
    (subnormal `x`), the Hessian is NaN for `0 < x <= 1e-155`, and the
    value is NaN for `x >= 9e307`, where `u` overflows
    (dev/famuni-rev2-out/guard.txt, hess.txt). The GEV and the Frechet
    cannot reach these (`1 + xi` and `nu - 1` are 0 or at least
    1.1e-16); frmtmb.eam's wiener guard reaches them only with a
    parameter below 1e-154."
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
