# Punch item 1, re-constructing reviewer findings 1 and 2 against the
# read-time link fields. Lane only.
#
#   Rscript dev/famlink-punch-view.R > dev/famlink-punch-view-log.txt
#
# Part A takes the reviewer's mutation paths from dev/famlink-rev-mutate.R
# and, for each, compares what `$` reports with the link the object's
# `links` holds, which is what a fit applies. Part B repeats the
# reviewer's fit construction (seed 20260916, n = 80). Part C reads the
# fits the BASE build saved (dev/famlink-rev-savedfit-base.rds, seed
# 20260916, n = 200) and runs `$link` and insight::model_info(). Part D
# counts calls to `$.frmtmb_family` during the reviewer's workload, with
# a positive control, instead of timing them.
ARM <- "lane"
source("dev/famlink-rev-common.R")
L <- frmtmb:::frmtmb_links
g <- function() brmsfamily("gaussian")

# the link each field should report, read straight from `links`
# (a missing link reads as NA, the value `$` gives for it)
truth <- function(f) {
  u <- unclass(f)
  nm <- function(lk) {
    v <- if (is.list(lk)) lk$name else lk
    if (is.null(v)) NA_character_ else v
  }
  list(link = nm(u$ord_link %||% u$links$mu),
       link_sigma = if ("sigma" %in% u$dpars) nm(u$links$sigma))
}
show <- function(label, expr) {
  v <- tryCatch(expr, error = function(e) e)
  if (inherits(v, "error")) {
    cat(sprintf("%-44s ERROR at construction: %s\n", label,
                substr(conditionMessage(v), 1, 60)))
    return(invisible())
  }
  got <- tryCatch(list(link = v$link, link_sigma = v$link_sigma),
                  error = function(e) e)
  if (inherits(got, "error")) {
    cat(sprintf("%-44s ERROR on read: %s\n", label,
                substr(conditionMessage(got), 1, 60)))
    return(invisible())
  }
  tr <- truth(v)
  ok <- identical(got$link, tr$link) &&
    identical(got$link_sigma, tr$link_sigma)
  cat(sprintf("%-44s %s  link=%s link_sigma=%s\n", label,
              if (ok) "agrees" else "DISAGREES", format(got$link),
              format(got$link_sigma)))
}

cat("-- A: mutation paths --\n")
show("fam$links$mu <- probit obj", { f <- g(); f$links$mu <- L$probit; f })
show("fam[['links']][['sigma']] <- 'softplus'",
     { f <- g(); f[["links"]][["sigma"]] <- "softplus"; f })
show("fam['links'] <- list(...)  ([<-)", {
  f <- g(); lk <- f[["links"]]; lk$mu <- L$log; lk$sigma <- L$softplus
  f["links"] <- list(lk); f })
show("fam[c('links','dpars')] <- ...", {
  f <- g(); lk <- f[["links"]]; lk$mu <- L$log
  f[c("links", "dpars")] <- list(lk, f[["dpars"]]); f })
show("fam[[idx of links]] <- ...", {
  f <- g(); lk <- f[["links"]]; lk$mu <- L$log
  f[[which(names(f) == "links")]] <- lk; f })
show("utils::modifyList(fam, list(links=))", {
  f <- g(); utils::modifyList(f, list(links = list(mu = L$log))) })
show("unclass, edit, class<-", {
  f <- g(); u <- unclass(f); u$links$mu <- L$log
  class(u) <- "frmtmb_family"; u })
show("structure(modified list, class=)", {
  f <- g(); u <- unclass(f); u$links$mu <- L$log
  structure(u, class = "frmtmb_family") })
show("fam$links <- NULL", { f <- g(); f$links <- NULL; f })
show("fam$ord_link on cumulative <- probit obj",
     { f <- cumulative(); f$ord_link <- L$probit; f })
show("fam['link'] <- 'log' ([<- store)", { f <- g(); f["link"] <- "log"; f })
show("fam$link <- 'log' (store)", { f <- g(); f$link <- "log"; f })
show("fam$link_sigma <- 'identity' (store)",
     { f <- g(); f$link_sigma <- "identity"; f })

cat("\n-- B: does family(fit) report the link the fit applied? --\n")
set.seed(20260916)
d <- data.frame(x = rnorm(80))
d$y <- exp(0.3 + 0.4 * d$x + rnorm(80, 0, 0.2))
f <- g()
lk <- f[["links"]]
lk$mu <- L$log
f["links"] <- list(lk)
fit <- frm(y ~ x, data = d, family = f)
cat("family(fit)$link:", family(fit)$link, " link applied by the fit:",
    fit$frame[["linpreds"]][[1]][["link"]][["name"]], "\n")

cat("\n-- C: fits saved by the base build --\n")
suppressMessages(library(insight))
fits <- readRDS("dev/famlink-rev-savedfit-base.rds")
for (nm in names(fits)) {
  fam <- family(fits[[nm]])
  lk <- tryCatch(fam$link, error = function(e) paste("ERROR:", conditionMessage(e)))
  mi <- tryCatch({
    m <- insight::model_info(fits[[nm]])
    paste0("ok, link_function = ", m$link_function)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-12s $link = %-10s model_info: %s\n", nm, lk, substr(mi, 1, 70)))
}

cat("\n-- D: calls to $.frmtmb_family, counted --\n")
set.seed(20260916)
d <- data.frame(x = rnorm(300), g = factor(rep(1:30, 10)))
d$y <- rpois(300, exp(0.5 + 0.3 * d$x + rnorm(30, 0, 0.3)[d$g]))
ns <- asNamespace("frmtmb")
orig <- get("$.frmtmb_family", ns)
cnt <- new.env()
cnt$n <- 0L
counting <- function(x, name) {
  cnt$n <- cnt$n + 1L
  orig(x, name)
}
unlockBinding("$.frmtmb_family", ns)
assign("$.frmtmb_family", counting, ns)
registerS3method("$", "frmtmb_family", counting, envir = ns)
tm <- function(label, expr) {
  cnt$n <- 0L
  force(expr)
  cat(sprintf("%-40s $ calls %7d\n", label, cnt$n))
}
tm("control: 3 direct reads (must be 3)",
   { ff <- negbinomial(); ff$lpdf; ff$links; ff$link })
tm("frm negbinomial + (1|g)",
   fit <- frm(y ~ x + (1 | g), data = d, family = negbinomial()))
tm("predict(newdata)", predict(fit, newdata = d))
tm("simulate(nsim = 50)", simulate(fit, nsim = 50, seed = 1))
tm("frm_bootstrap(nsim = 20)", frm_bootstrap(fit, nsim = 20, seed = 1))
tm("print(family(fit))", capture.output(print(family(fit))))
assign("$.frmtmb_family", orig, ns)
registerS3method("$", "frmtmb_family", orig, envir = ns)
