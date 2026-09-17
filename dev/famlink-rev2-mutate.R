## Recheck (rounds 1 and 1b), priority 2: read-time link fields.
## The round 0 mutation paths (dev/famlink-rev-mutate.R) plus more, each
## judged against an INDEPENDENT reading of the object: the `name`,
## `linkfun` and `linkinv` of the link objects in `links` (or
## `ord_link`), read with .subset2 so the method under test does not
## judge itself. Deterministic; the fit at the end uses seed 20260916,
## n = 80. Usage: Rscript dev/famlink-rev2-mutate.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
L <- frmtmb:::frmtmb_links
truth_name <- function(lk) {
  nm <- if (is.list(lk)) lk[["name"]] else lk
  if (is.character(nm) && length(nm) == 1L) nm else NA_character_
}
dollar <- function(f, nm) eval(call("$", f, nm))
judge <- function(f) {
  lk <- .subset2(f, "links"); ol <- .subset2(f, "ord_link")
  dp <- .subset2(f, "dpars")
  joint <- identical(.subset2(f, "type"), "categorical") ||
    identical(.subset2(f, "family"), "multinomial")
  main <- if (!is.null(ol)) ol else if (joint) "logit" else
    if ("mu" %in% dp) lk[["mu"]] else "identity"
  bad <- character()
  if (!identical(f$link, truth_name(main))) bad <- c(bad, "link")
  mobj <- if (is.list(main)) main else
    tryCatch(frmtmb:::get_link(main), error = function(e) NULL)
  if (!identical(f$linkinv, mobj[["linkinv"]])) bad <- c(bad, "linkinv")
  if (!identical(f$linkfun, mobj[["linkfun"]])) bad <- c(bad, "linkfun")
  sec <- setdiff(dp, "mu")
  if (joint) sec <- setdiff(sec, .subset2(f, "primary_dpars"))
  for (d in sec) {
    if (!identical(dollar(f, paste0("link_", d)), truth_name(lk[[d]]))) {
      bad <- c(bad, paste0("link_", d))
    }
  }
  if (length(bad)) paste("DISAGREES:", paste(bad, collapse = ",")) else
    paste0("agrees (", 3L + length(sec), " fields)")
}
show <- function(label, expr) {
  v <- tryCatch(expr, error = function(e) e)
  if (inherits(v, "error")) {
    cat(sprintf("%-42s ERROR building: %s\n", label,
                substr(conditionMessage(v), 1, 70)))
    return(invisible())
  }
  j <- tryCatch(judge(v), error = function(e) {
    paste("ERROR on read:", substr(conditionMessage(e), 1, 60))
  })
  cat(sprintf("%-42s %s  link=%s\n", label, j,
              tryCatch(format(v$link), error = function(e) "<error>")))
}
g <- function() brmsfamily("gaussian")
show("fam$links$mu <- probit obj", { f <- g(); f$links$mu <- L$probit; f })
show("fam[['links']][['sigma']] <- 'softplus'", { f <- g(); f[["links"]][["sigma"]] <- "softplus"; f })
show("fam['links'] <- list(...)", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; lk$sigma <- L$softplus; f["links"] <- list(lk); f })
show("fam[c('links','dpars')] <- ...", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f[c("links", "dpars")] <- list(lk, f[["dpars"]]); f })
show("fam[[idx of links]] <- ...", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f[[which(names(f) == "links")]] <- lk; f })
show("utils::modifyList(fam, list(links=))", { f <- g(); utils::modifyList(f, list(links = list(mu = L$log))) })
show("unclass, edit, class<-", { f <- g(); u <- unclass(f); u$links$mu <- L$log; class(u) <- "frmtmb_family"; u })
show("structure(modified, class=)", { f <- g(); u <- unclass(f); u$links$mu <- L$log; structure(u, class = "frmtmb_family") })
show("fam$dpars <- c(dpars, 'extra')", { f <- g(); f$dpars <- c(f$dpars, "extra"); f })
show("fam$links <- NULL", { f <- g(); f$links <- NULL; f })
show("fam$ord_link on cumulative <- probit", { f <- cumulative(); f$ord_link <- L$probit; f })
show("sratio('cloglog')", sratio("cloglog"))
show("mixture(gaussian, student)", mixture(gaussian, student))
show("mixture, component link edited", { f <- mixture(gaussian, gaussian); f$links$mu2 <- L$log; f })
show("categorical(levels = a, b, c)", categorical(levels = c("a", "b", "c")))
show("multinomial(3)", multinomial(3))
show("custom family, string links", frmtmb_family("cf", dpars = c("mu", "phi"), links = list(mu = "log", phi = "softplus"), lpdf = function(y, dpars, aterms) 0))
show("custom family, unknown link string", frmtmb_family("cf", dpars = "mu", links = list(mu = "nosuch"), lpdf = function(y, dpars, aterms) 0))
show("fam['link'] <- 'log' (store)", { f <- g(); f["link"] <- "log"; f })
show("fam$link_sigma <- 'identity' (store)", { f <- g(); f$link_sigma <- "identity"; f })
cat("\n-- reads that do not go through `$` --\n")
f <- g()
cat("f[['link']]:", format(f[["link"]]), "| unclass(f)$link:",
    format(unclass(f)$link), "| getElement(f, 'link'):",
    format(getElement(f, "link")), "| 'link' %in% names(f):",
    "link" %in% names(f), "\n")
h <- g(); h$links <- NULL
cat("links removed: $link", format(h$link), " $linkinv is NULL:",
    is.null(h$linkinv), " $link_sigma", format(h$link_sigma), "\n")
cat("\n-- through a fit --\n")
set.seed(20260916)
d <- data.frame(x = rnorm(80)); d$y <- exp(0.3 + 0.4 * d$x + rnorm(80, 0, 0.2))
f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f["links"] <- list(lk)
fit <- frm(y ~ x, data = d, family = f)
cat("family(fit)$link:", family(fit)$link, " link applied by the fit:",
    fit$frame$linpreds[[1]]$link$name, "\n")
