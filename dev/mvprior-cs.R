# Lane wt-mvprior, punch round 2: class "b" and cs() on the sequential
# ordinal families. brms 2.23.0's stancode() lines and default_prior()
# rows, and frmtmb's resolved parameters per arm, on the same data.
#   MVPRIOR_WHO=brms Rscript dev/mvprior-cs.R > dev/mvprior-log/cs-brms.txt
#   MVPRIOR_ARM=base|lane Rscript dev/mvprior-cs.R > .../cs-<arm>.txt
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
who <- Sys.getenv("MVPRIOR_WHO", mvprior_arm)
set.seed(4417)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
d$ord <- factor(cut(d$x + 0.5 * d$z + rlogis(n), c(-Inf, -1, 0, 1.2, Inf),
                    labels = FALSE), ordered = TRUE)
forms <- list(ord ~ cs(x), ord ~ z + cs(x), ord ~ z + cs(x) + cs(w))
specs <- list(list(class = "b", coef = ""), list(class = "b", coef = "x"),
              list(class = "b", coef = "z"), list(class = "b", coef = "w"))
flat <- function(x) gsub("[[:space:]]+", " ", x)
if (identical(who, "brms")) {
  suppressMessages(library(brms))
  cat("brms", format(packageVersion("brms")), "\n")
} else {
  suppressMessages(library(frmtmb))
  cat("arm", mvprior_arm, "frmtmb from", find.package("frmtmb"), "\n")
  ns <- asNamespace("frmtmb")
}
for (famname in c("acat", "sratio", "cratio")) {
  for (f in forms) {
    cat("\n==", famname, deparse(f), "==\n")
    for (s in specs) {
      lab <- sprintf("b coef=%-2s", s$coef)
      if (identical(who, "brms")) {
        fam <- get(famname, asNamespace("brms"))()
        r <- tryCatch({
          sc <- strsplit(stancode(f, data = d, family = fam,
                                  prior = brms::set_prior("normal(0, 0.123)",
                                                          class = s$class,
                                                          coef = s$coef)),
                         "\n")[[1]]
          paste(trimws(grep("0.123", sc, value = TRUE)), collapse = " | ")
        }, error = function(e) paste("REFUSE:", flat(conditionMessage(e))))
      } else {
        fam <- get(famname, asNamespace("frmtmb"))()
        des <- ns$prior_design(f, d, fam, list())
        pt <- des$frame[["par_template"]]
        r <- tryCatch({
          rr <- ns$resolve_priorlist(des, set_prior("normal(0, 0.123)",
                                                    class = s$class,
                                                    coef = s$coef))
          nm <- vapply(rr$entries, function(e) {
            paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx],
                  collapse = "+")
          }, "")
          if (length(nm)) paste(nm, collapse = " ") else "(no density)"
        }, error = function(e) paste("REFUSE:", flat(conditionMessage(e))))
      }
      cat(lab, r, "\n")
    }
    tab <- if (identical(who, "brms")) {
      default_prior(f, data = d, family = fam)
    } else {
      default_prior(f, data = d, family = fam)
    }
    tab <- as.data.frame(tab)
    tab <- tab[tab$class == "b", ]
    cat("b rows:", paste0("[", tab$coef, "]", collapse = " "), "\n")
  }
}
