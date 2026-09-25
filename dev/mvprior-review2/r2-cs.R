# cs() on an ordinal family: class "b" with no ordinary slope.
who <- Sys.getenv("R2_WHO", "lane")
Sys.setenv(R2_ARM = if (who == "base") "base" else "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
if (who == "brms") {
  suppressPackageStartupMessages(library(brms))
  sc <- stancode(ord ~ cs(x), data = d, family = acat(),
                 prior = set_prior("normal(0, 0.123)", class = "b"))
  cat(grep("0.123|bcs", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  print(default_prior(ord ~ cs(x), data = d, family = acat()))
  sc2 <- stancode(ord ~ z + cs(x), data = d, family = acat(),
                  prior = set_prior("normal(0, 0.123)", class = "b"))
  cat("--- ord ~ z + cs(x)\n")
  cat(grep("0.123", strsplit(sc2, "\n")[[1]], value = TRUE), sep = "\n")
} else {
  ns <- asNamespace("frmtmb")
  for (f in list(ord ~ cs(x), ord ~ z + cs(x))) {
    des <- ns$prior_design(f, d, acat(), list())
    pt <- des$frame[["par_template"]]
    r <- tryCatch({
      rr <- ns$resolve_priorlist(des, set_prior("normal(0, 0.123)", class = "b"))
      paste(vapply(rr$entries, function(e) paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+"), ""), collapse = " ")
    }, error = function(e) paste("REFUSE:", conditionMessage(e)))
    cat(deparse(f), "->", if (nzchar(r)) r else "(no density)", "\n")
    cat("  par names:", unlist(lapply(names(pt), function(c) ns$par_template_names(pt[[c]], c))), "\n")
    print(as.data.frame(default_prior(f, data = d, family = acat()))[, c("class", "coef", "dpar")])
  }
}
