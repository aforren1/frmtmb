# Short check, item 7: resolve the refreshed trial merge. Each textual hunk keeps both
# lanes' intent; S2, S4 and S5 of dev/predread-findings.md section 8 are applied as exact
# edits, each asserted to match once. S3 needs no edit: latres now keeps dim(v).
M <- "C:/Users/adf44/source/r/rev-predread-src/merge"
res <- function(f, pick) {
  p <- file.path(M, f)
  x <- readLines(p, warn = FALSE)
  out <- character()
  i <- 1L
  n <- 0L
  while (i <= length(x)) {
    if (startsWith(x[i], "<<<<<<< ")) {
      rest <- (i + 1L):length(x)
      j1 <- rest[startsWith(x[rest], "||||||| ")][1L]
      j2 <- rest[startsWith(x[rest], "=======")][1L]
      j3 <- rest[startsWith(x[rest], ">>>>>>> ")][1L]
      ours <- if (j1 > i + 1L) x[(i + 1L):(j1 - 1L)] else character()
      theirs <- if (j3 > j2 + 1L) x[(j2 + 1L):(j3 - 1L)] else character()
      n <- n + 1L
      out <- c(out, pick(ours, theirs, n))
      i <- j3 + 1L
    } else {
      out <- c(out, x[i])
      i <- i + 1L
    }
  }
  writeLines(out, p)
  cat(f, "hunks:", n, "\n")
}
sub1 <- function(f, old, new) {
  p <- file.path(M, f)
  x <- paste(readLines(p, warn = FALSE), collapse = "\n")
  parts <- strsplit(x, old, fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("pattern matched ", length(parts) - 1L, " times in ", f, ": ", old)
  writeLines(paste0(parts[1L], new, parts[2L]), p, useBytes = TRUE)
  cat(f, "edited\n")
}
# fitted_point_se(): latres's early returns, on no me() row, then predread's test
res("R/predict.R", function(ours, theirs, n) {
  t2 <- theirs[-length(theirs)]
  t2 <- sub("if (!is.null(newdata) && !is.matrix(est) && scale == \"response\" &&",
            "if (is.null(mm) && !is.null(newdata) && !is.matrix(est) && scale == \"response\" &&",
            t2, fixed = TRUE)
  c(t2, ours)
})
# get_predict(): one mean, S4
res("R/interop.R", function(ours, theirs, n) {
  c("  # the expected response over a new row's me() value and latent residual",
    "  if (type == \"response\") {",
    "    mm <- fitted_me_moments(model, newdata, dots[[\"re_formula\"]],",
    "                            dots[[\"resp\"]], dpar,",
    "                            isTRUE(dots[[\"allow_new_levels\"]]))",
    "    if (!is.null(mm)) {",
    "      attributes(mm$mean) <- attributes(p)",
    "      p <- mm$mean",
    "    } else {",
    "      p <- latres_fitted_mean(model, newdata, dots[[\"re_formula\"]],",
    "                              dots[[\"resp\"]], dpar,",
    "                              isTRUE(dots[[\"allow_new_levels\"]]), p)",
    "    }",
    "  }")
})
res("extensions/frmtmb.sample/R/methods-draws.R", function(ours, theirs, n) {
  c(theirs, ours[-1L])
})
# S5 and the predict row: nlfix's nl row; latres's first sentence on predread's row
res("R/compat.R", function(ours, theirs, n) {
  if (n == 1L) return(theirs)
  sub("\"Under cov = TRUE predict() draws",
      "\"On a family with no residual y - mu, predict() reads the latent's conditional mode on fitted rows and draws new latents on new rows under cov = TRUE. On gaussian and student, under cov = TRUE predict() draws",
      ours, fixed = TRUE)
})
for (f in c("NEWS.md", "extensions/frmtmb.sample/NEWS.md")) {
  res(f, function(ours, theirs, n) c(ours, theirs))
}
# S2
sub1("R/predict.R",
     "out <- if (type == \"response\") lp[[\"link\"]]$linkinv(eta) else eta\n      return(if (is.null(newdata)) napred(object, out) else out)",
     "out <- if (type == \"response\") {\n        dpar_response_link(rspec, lp, dpar)$linkinv(eta)\n      } else {\n        eta\n      }\n      return(if (is.null(newdata)) napred(object, out) else out)")
# S4: one response_mean_gh(), latres's, which has the tail guard
p <- file.path(M, "R/me.R")
x <- readLines(p, warn = FALSE)
st <- which(startsWith(x, "response_mean_gh <- function"))
stopifnot(length(st) == 1L)
en <- st + which(x[(st + 1L):length(x)] == "}")[1L]
x <- x[-(st:en)]
writeLines(x, p)
cat("R/me.R: response_mean_gh removed, lines", st, "to", en, "\n")
# S4: the latent residual's variance joins the node when the fit also has me()
sub1("R/predict.R",
     "  if (is.null(loads)) return(NULL)\n  lp <- object$frame[[\"linpreds\"]][[loads$key]]",
     paste0("  if (is.null(loads)) return(NULL)\n",
            "  bk <- frame_latres(object$frame)[[rn]]\n",
            "  if (!is.null(bk) && !isTRUE(object$frame[[\"latres_off\"]]) &&\n",
            "        identical(loads$key, linpred_key(rn, \"mu\"))) {\n",
            "    sl <- latres_marginal_sd(object$estimates[[\"theta\"]][bk[[\"theta_idx\"]]],\n",
            "                             bk[[\"aux\"]])\n",
            "    s_old <- sqrt(loads$s2)\n",
            "    loads$s2 <- loads$s2 + sl^2\n",
            "    s_new <- sqrt(loads$s2)\n",
            "    for (k in names(loads$by_key)) {\n",
            "      loads$by_key[[k]] <- ifelse(s_new > 0, loads$by_key[[k]] * s_old / s_new, 0)\n",
            "    }\n",
            "    loads$by_key[[loads$key]] <- s_new\n",
            "  }\n",
            "  lp <- object$frame[[\"linpreds\"]][[loads$key]]"))
sub1("R/predict.R",
     "  if (!is.null(newdata) && scale == \"response\" &&\n        (is.null(dpar) || identical(dpar, \"mu\"))) {\n    out <- latres_fitted_mean(",
     "  if (is.null(mm) && !is.null(newdata) && scale == \"response\" &&\n        (is.null(dpar) || identical(dpar, \"mu\"))) {\n    out <- latres_fitted_mean(")
left <- system2("grep", c("-l", shQuote("^<<<<<<< "), "-r", shQuote(file.path(M, "R")),
                          shQuote(file.path(M, "extensions/frmtmb.sample/R"))),
                stdout = TRUE)
cat("files with markers left:", length(left), "\n")
