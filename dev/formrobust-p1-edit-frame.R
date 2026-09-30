# Punch round 1, minors m5, m6, m7 (cens y2), m8: record of the edit of
# R/frame.R.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/frame.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches for ", substr(old, 1, 70))
  x <<- sub(old, new, x, fixed = TRUE)
}
# m6: a bare name reads a scalar constant of the formula environment as
# the constant it is, the way trials(k + 0) and trials(10) already did
rep1("      if (is.name(a)) {
        add_part(a)
        next
      }",
"      if (is.name(a)) {
        # a scalar of the formula environment is the constant it names,
        # as the literal trials(10) is, so it stays out of the frame
        for (v in expr_frame_vars(a, data, resp$formula_env)) {
          add_part(as.name(v))
        }
        next
      }")
# m5: an NA the expression produces is refused, as brms refuses it
rep1("          if (length(v) == 1L) v <- rep(v, n)
          if (!is.matrix(v) && length(v) != n) {",
"          if (length(v) == 1L) v <- rep(v, n)
          if (anyNA(v)) {
            # the variables are complete by here, so the expression made
            # the NA; brms refuses it at standata(), and a row dropped
            # for it would change the sample in silence
            frm_stop(\"Addition term \", aterm_label(nm_at, a), \" of \",
                     \"response '\", resp$resp_name, \"' is NA on \",
                     sum(is.na(v)), \" of \", n, \" rows, where its \",
                     \"variables are not. Give every row a value\",
                     call. = FALSE)
          }
          if (!is.matrix(v) && length(v) != n) {")
# m7: a single-value interval bound is refused in brms's sense rather
# than read as NA on the censored rows
rep1("      v <- as.numeric(eval(resp$aterms[[\"cens_y2\"]], data, resp$formula_env))
      if (!is.null(attr(mf, \"na.action\"))) {",
"      v <- as.numeric(eval(resp$aterms[[\"cens_y2\"]], data, resp$formula_env))
      if (length(v) != NROW(data)) {
        # brms: \"Argument 'y2' needs to have length equal to the number
        # of data rows\"; it is not recycled, unlike the censoring code
        frm_stop(\"cens(\", deparse1(resp$aterms[[\"cens\"]]), \", \",
                 deparse1(resp$aterms[[\"cens_y2\"]]), \"): the interval \",
                 \"upper bound has \", length(v), \" value(s) where the data \",
                 \"have \", NROW(data), \" rows. It takes one value per row, \",
                 \"as brms requires\", call. = FALSE)
      }
      if (!is.null(attr(mf, \"na.action\"))) {")
# m8: a name that R finds as a function is not a variable
rep1("    x <- tryCatch(get(v, envir = env), error = function(e) NULL)
    # an object that is not found goes to the frame so that the frame's
    # own check names it; a constant such as `k` in weights(wt * k) is
    # read when the term is evaluated, since it has no row to drop",
"    x <- tryCatch(get(v, envir = env), error = function(e) NULL)
    if (is.function(x)) {
      # weights(t * 2) without a column t found base::t(), and the
      # arithmetic then failed without naming the variable
      frm_stop(\"The model uses `\", v, \"`, which is not a column of \",
               \"`data`; R finds only the function \", v, \"() from the \",
               \"formula. Add the column to `data` or correct the name\",
               call. = FALSE)
    }
    # an object that is not found goes to the frame so that the frame's
    # own check names it; a constant such as `k` in weights(wt * k) is
    # read when the term is evaluated, since it has no row to drop")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
