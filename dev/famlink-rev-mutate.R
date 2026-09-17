## Reviewer check for lane wt-famlink, priority 1: every way a user or
## an extension can change `links` after construction, and whether
## `link` / `link_<dpar>` follow. Deterministic, no data.
## Usage: Rscript dev/famlink-rev-mutate.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
L <- frmtmb:::frmtmb_links
stale <- function(f) {
  u <- unclass(f)
  fresh <- unclass(frmtmb:::family_link_fields(f))
  nm <- union(grep("^link", names(u), value = TRUE), grep("^link", names(fresh), value = TRUE))
  nm <- setdiff(nm, c("links", "linkfun", "linkinv"))
  bad <- nm[!vapply(nm, function(k) identical(u[[k]], fresh[[k]]), TRUE)]
  li <- tryCatch(!isTRUE(all.equal(u[["linkinv"]](0.3), fresh[["linkinv"]](0.3))),
                 error = function(e) NA)
  paste0(if (length(bad)) paste(bad, collapse = ",") else "fresh",
         if (isTRUE(li)) " +linkinv stale" else "")
}
show <- function(label, expr) {
  v <- tryCatch(expr, error = function(e) e)
  cat(sprintf("%-46s %s\n", label, if (inherits(v, "error"))
    paste("ERROR:", substr(conditionMessage(v), 1, 70)) else
    paste0(stale(v), "  link=", unclass(v)$link,
           "  link_sigma=", format(unclass(v)$link_sigma))))
}
g <- function() brmsfamily("gaussian")
show("fam$links$mu <- probit obj", { f <- g(); f$links$mu <- L$probit; f })
show("fam[['links']][['sigma']] <- 'softplus'", { f <- g(); f[["links"]][["sigma"]] <- "softplus"; f })
show("fam['links'] <- list(...)  ([<-)", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; lk$sigma <- L$softplus; f["links"] <- list(lk); f })
show("fam[c('links','dpars')] <- ...", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f[c("links", "dpars")] <- list(lk, f[["dpars"]]); f })
show("fam[[idx of links]] <- ...", { f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f[[which(names(f) == "links")]] <- lk; f })
show("utils::modifyList(fam, list(links=))", { f <- g(); utils::modifyList(f, list(links = list(mu = L$log))) })
show("unclass, edit, class<-", { f <- g(); u <- unclass(f); u$links$mu <- L$log; class(u) <- "frmtmb_family"; u })
show("structure(modified list, class=)", { f <- g(); u <- unclass(f); u$links$mu <- L$log; structure(u, class = "frmtmb_family") })
show("within(fam, links$mu <- ...)", { f <- g(); within(unclass(f), links$mu <- L$log) |> structure(class = "frmtmb_family") })
show("attr<- names (rename a field)", { f <- g(); names(f)[names(f) == "links"] <- "links"; f })
show("fam$dpars <- c(dpars, 'extra')", { f <- g(); f$dpars <- c(f$dpars, "extra"); f })
show("fam$links <- NULL", { f <- g(); f$links <- NULL; f })
show("fam$ord_link on cumulative <- probit obj", { f <- cumulative(); f$ord_link <- L$probit; f })
show("mixture comps edited: fam$mix$... ", { f <- mixture(gaussian, gaussian); f$mix$K <- 2; f })
show("fam$link <- 'log' (refused?)", { f <- g(); f$link <- "log"; f })
show("fam[['link']] <- 'log' (refused?)", { f <- g(); f[["link"]] <- "log"; f })
show("fam['link'] <- 'log' ([<- bypass)", { f <- g(); f["link"] <- "log"; f })
show("fam$linkinv <- exp (refused?)", { f <- g(); f$linkinv <- exp; f })
show("fam$link_sigma <- 'identity' (refused?)", { f <- g(); f$link_sigma <- "identity"; f })
show("fam$link_new <- 'x' on absent dpar", { f <- g(); f$link_new <- "x"; f })
cat("\n-- does a stale `link` after [<- change the fit? --\n")
set.seed(20260916)
d <- data.frame(x = rnorm(80)); d$y <- exp(0.3 + 0.4 * d$x + rnorm(80, 0, 0.2))
f <- g(); lk <- f[["links"]]; lk$mu <- L$log; f["links"] <- list(lk)
fit <- frm(y ~ x, data = d, family = f)
cat("family(fit)$link:", family(fit)$link, " links$mu used by the fit:",
    fit$frame$linpreds[[1]]$link$name, "\n")
