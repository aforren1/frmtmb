# Reviewer check (lane ceplot): what in a fit changes when
# conditional_effects() runs on it? Every environment reachable from the
# fit is listed with a hash of its contents before and after, for a
# relabel call (unseen g:h) and for a control call (observed g:h), on
# the fit (boot) and on draws.
#   Rscript dev/ceplot-rev-leak.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
envs <- function(x, path = "fit", acc = new.env()) {
  if (is.environment(x)) {
    if (identical(x, globalenv()) || identical(x, emptyenv()) ||
        isNamespace(x) || identical(x, baseenv())) return(acc)
    key <- format(x)
    if (!is.null(acc[[key]])) return(acc)
    acc[[key]] <- list(path = path, env = x)
    for (n in ls(x, all.names = TRUE)) {
      v <- tryCatch(get(n, envir = x), error = function(e) NULL)
      envs(v, paste0(path, "$", n), acc)
    }
    return(acc)
  }
  if (is.function(x)) return(envs(environment(x), paste0(path, "<fenv>"), acc))
  if (inherits(x, "formula")) {
    return(envs(environment(x), paste0(path, "<fmlenv>"), acc))
  }
  if (is.list(x)) {
    for (i in seq_along(x)) {
      nm <- names(x)[i]
      if (is.null(nm) || !nzchar(nm)) nm <- paste0("[[", i, "]]")
      envs(x[[i]], paste0(path, "$", nm), acc)
    }
  }
  a <- attributes(x)
  if (!is.null(a) && !is.list(x)) {
    for (n in names(a)) envs(a[[n]], paste0(path, "@", n), acc)
  }
  acc
}
snap <- function(fit) {
  acc <- envs(fit)
  out <- lapply(ls(acc), function(k) {
    e <- acc[[k]]$env
    list(path = acc[[k]]$path, names = sort(ls(e, all.names = TRUE)),
         hash = digest_env(e))
  })
  names(out) <- vapply(out, `[[`, "", "path")
  out
}
digest_env <- function(e) {
  vals <- lapply(sort(ls(e, all.names = TRUE)), function(n) {
    v <- tryCatch(get(n, envir = e), error = function(err) "<unbound>")
    if (is.environment(v)) "<env>" else if (is.function(v)) {
      deparse(v)
    } else v
  })
  paste(as.character(tools::md5sum(textConnection_to_file(vals))))
}
textConnection_to_file <- function(v) {
  f <- tempfile()
  saveRDS(v, f, compress = FALSE)
  f
}
compare <- function(a, b, what) {
  ch <- character(0)
  for (p in union(names(a), names(b))) {
    if (is.null(a[[p]]) || is.null(b[[p]])) {
      ch <- c(ch, paste0(p, " (appeared/vanished)"))
    } else if (!identical(a[[p]]$hash, b[[p]]$hash)) {
      nn <- setdiff(b[[p]]$names, a[[p]]$names)
      ch <- c(ch, paste0(p, " contents changed; new names: ",
                         paste(nn, collapse = ",")))
    }
  }
  cat(what, ":", if (length(ch)) "" else "no environment changed", "\n")
  for (c1 in ch) cat("   ", c1, "\n")
}
S <- shape_fits()
fit <- S$A$fit
fr0 <- fit$frame
es0 <- fit$estimates
ce1 <- function(o, cond, ...) {
  suppressWarnings(conditional_effects(o, "x", resolution = 3,
                                       re_formula = NULL, conditions = cond,
                                       ...))
}
ds <- hand(fit, 50)
invisible(ce1(fit, list(g = "1", h = "2")))  # warm the caches first
s0 <- snap(fit)
invisible(ce1(fit, list(g = "1", h = "2"), band = "boot", boot = 10, seed = 1))
s1 <- snap(fit)
compare(s0, s1, "control (observed g:h), boot")
r <- tryCatch(ce1(fit, list(g = "1", h = "1"), band = "boot", boot = 10,
                  seed = 1), error = function(e) conditionMessage(e))
s2 <- snap(fit)
compare(s1, s2, "unseen g:h, boot")
invisible(ce1(ds, list(g = "1", h = "2"), seed = 1))
s3 <- snap(fit)
compare(s2, s3, "control on draws")
r <- tryCatch(ce1(ds, list(g = "1", h = "1"), seed = 1),
              error = function(e) conditionMessage(e))
s4 <- snap(fit)
compare(s3, s4, "unseen g:h on draws")
cat("frame identical to the value before:", identical(fr0, fit$frame), "\n")
cat("estimates identical:", identical(es0, fit$estimates), "\n")
lv <- lapply(fit$frame$re_blocks, `[[`, "levels")
cat("block levels contain '1:1':", any(vapply(lv, function(l) "1:1" %in% l,
                                                NA)), "\n")
# the fit's own predictions afterwards equal a fresh fit's
fresh <- shape_fits()$A$fit
nd <- data.frame(x = c(0, 1), g = factor(c("1", "2"), levels = 1:6),
                 h = factor(c("2", "1"), levels = 1:5))
cat("fitted() after the calls identical to a fresh fit's:",
    identical(fitted(fit, newdata = nd), fitted(fresh, newdata = nd)), "\n")
cat("wald ce identical to a fresh fit's:",
    identical(ce1(fit, list(g = "1", h = "1"))$x$estimate__,
              ce1(fresh, list(g = "1", h = "1"))$x$estimate__), "\n")
