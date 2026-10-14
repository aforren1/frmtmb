# Release reviewer of 0.71.0: the error of each postfix stress model
# (dev/postfix-fix2-nlnew.R, "stress: sinusoid") on rellib-r9, checked
# for nl_start_message()'s text and for any scan note.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
f1 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
lead <- "The nonlinear fit failed from its default starting values ("
tab <- NULL
for (P0 in c(1.3, 3, 4, 7, 9)) for (s in 1:6) {
  set.seed(1000 + s)
  xs <- seq(0, 15, length.out = 120)
  ds <- data.frame(x = xs, y = 2 * sin(6.283185307179586 * xs / P0) +
                     rnorm(120, 0, 0.6))
  e <- tryCatch({ suppressWarnings(frm(f1, data = ds)); NULL },
                error = function(e) e)
  m <- if (is.null(e)) "" else conditionMessage(e)
  tab <- rbind(tab, data.frame(P = P0, seed = 1000 + s, error = !is.null(e),
    nl_start_message = startsWith(m, lead),
    names_A_P = grepl("`A`, `P`", m, fixed = TRUE),
    scan = grepl("scan", m, ignore.case = TRUE),
    class = if (is.null(e)) "" else paste(class(e), collapse = ",")))
}
print(table(error = tab$error, nl_msg = tab$nl_start_message,
            names = tab$names_A_P, scan = tab$scan))
cat("classes:", unique(tab$class), "\n")
cat("one message:\n", m, "\n")
