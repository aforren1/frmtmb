p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
old <- "  # the objective is -logLik; a positive unit rescales, not the sign
  -g[j] / lp[[\"link\"]]$mu_eta(eta)"
new <- "  # the objective is -logLik, judged in the natural units the optimizer
  # used, as grad_verdict() judges it
  u <- fit$par_units %||% rep(1, length(p))
  -g[j] * u[j] / lp[[\"link\"]]$mu_eta(eta)"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
