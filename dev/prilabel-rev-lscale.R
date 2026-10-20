# Review of lane prilabel: a class-wide user prior over a default that
# frmtmb.sample writes per coefficient (lscale of a gp() term). Does
# the user's prior still land, and what does brms 2.23.0 do?
#   arg 1: "lane" or "base". Data seed 51.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
sink(file.path(wt, paste0("dev/prilabel-rev-log/lscale-", arm, ".txt")),
     split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "from",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(51)
n <- 80
d <- data.frame(x = runif(n), g = factor(sample(1:8, n, TRUE)))
d$y <- sin(4 * d$x) + rnorm(8, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)
uf <- frm(bf(y ~ gp(x) + (1 | g)), data = d, dry_run = "objective")
show <- function(pr) {
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, pr))$ri
  for (e in ri$entries) {
    cat(sprintf("  %-6s idx %-3s %-10s %s\n", e$comp, paste(e$idx, collapse = ","),
                e$dist$kind, paste(unlist(e$dist[setdiff(names(e$dist),
                                                        "kind")]),
                                   collapse = " ")))
  }
}
cat("\n-- defaults only\n")
show(NULL)
cat("\n-- user set_prior(\"inv_gamma(2, 1)\", class = \"lscale\")\n")
show(set_prior("inv_gamma(2, 1)", class = "lscale"))
cat("\n-- user set_prior(\"normal(0, 1)\", class = \"sd\", group = \"g\")\n")
show(set_prior("normal(0, 1)", class = "sd", group = "g"))
if (arm == "lane") {
  suppressPackageStartupMessages(library(brms))
  cat("\n-- brms: class-wide lscale prior beside its per-coef default\n")
  sc <- brms::stancode(y ~ gp(x) + (1 | g), data = d,
                       prior = brms::set_prior("inv_gamma(2, 1)",
                                               class = "lscale"))
  print(grep("lprior.*lscale", strsplit(sc, "\n")[[1]], value = TRUE))
  p <- as.data.frame(brms::default_prior(y ~ gp(x) + (1 | g), data = d))
  print(p[p$class == "lscale", c("prior", "class", "coef")])
}
sink()
