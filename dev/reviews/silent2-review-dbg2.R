# Review of lane silent2: what re_term_contrasts() sees for (1 + factor(cyl) | g).
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib", "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(7)
dc <- data.frame(g = factor(rep(1:10, each = 9)), x = rnorm(90),
                 cyl = sample(c(4, 6, 8), 90, TRUE))
dc$y <- rnorm(10)[dc$g] + rnorm(90)
trace("re_term_contrasts", where = asNamespace("frmtmb"), print = FALSE,
      tracer = quote({
        cat("re_term_contrasts lhs:", deparse(lhs), " data cols:",
            paste(names(data), collapse = ","), "\n")
        tt0 <- stats::terms(stats::as.formula(call("~", lhs), env = env))
        r <- tryCatch(stats::model.frame(tt0, data, na.action = stats::na.pass),
                      error = function(e) conditionMessage(e))
        cat("  model.frame:", if (is.character(r)) r else paste(names(r), collapse = ","), "\n")
      }))
fc <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
  suppressWarnings(suppressMessages(frm(y ~ (1 + factor(cyl) | g), data = dc))))
dc$fcyl <- factor(dc$cyl)
fc2 <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
  suppressWarnings(suppressMessages(frm(y ~ (1 + fcyl | g), data = dc))))
