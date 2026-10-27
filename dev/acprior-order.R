# Row order of frmtmb's default_prior() (after) against brms 2.23.0's,
# case by case, read from the two logs. Log: dev/acprior-out/order-vs-brms.txt
out <- "C:/Users/adf44/source/r/frmtmb-wt-acprior/dev/acprior-out"
rd <- function(f) {
  l <- readLines(file.path(out, f))
  res <- list()
  cur <- NULL
  for (x in l) {
    if (grepl("user priors|two slots", x)) {
      cur <- NULL
      next
    }
    if (grepl("^=====", x)) {
      cur <- sub("^===== ([^ ]+).*", "\\1", x)
      res[[cur]] <- character(0)
      next
    }
    if (is.null(cur) || !grepl("^ *[0-9]+ ", x)) next
    tk <- strsplit(trimws(x), " +")[[1]]
    cls <- intersect(tk, c("ar", "ma", "b", "Intercept", "sigma", "cosy",
                           "cortime", "theta"))[1]
    i <- match(cls, tk)
    co <- if (identical(cls, "b") && length(tk) > i && tk[i + 1] == "x") {
      "x"
    } else ""
    res[[cur]] <- c(res[[cur]], trimws(paste(cls, co)))
  }
  res
}
b <- rd("brms-defaults.txt")
a <- rd("frm-table-after.txt")
for (k in names(b)) {
  cat(sprintf("%-8s same=%s brms: %s | frmtmb: %s\n", k,
              identical(b[[k]], a[[k]]), paste(b[[k]], collapse = ", "),
              paste(a[[k]], collapse = ", ")))
}
