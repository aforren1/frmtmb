# Reviewer 2 probe. R2_WHO = brms | lane | base.
# For each case: default_prior() rows; for each spec: ACCEPT with the
# rows validate_prior() marks source = "user" (and, for frmtmb, the
# internal parameters the resolver reaches), or REFUSE with the message.
who <- Sys.getenv("R2_WHO", "lane")
Sys.setenv(R2_ARM = if (who == "base") "base" else "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
if (who == "brms") suppressPackageStartupMessages(library(brms))
ns_name <- if (who == "brms") "brms" else "frmtmb"
ns <- asNamespace(ns_name)
d <- r2_data()
flat <- function(x) gsub("[[:space:]]+", " ", x)
lab <- function(s) {
  f <- c(class = s$class, coef = s$coef, group = s$group, resp = s$resp,
         dpar = s$dpar, nlpar = s$nlpar)
  f <- f[nzchar(f)]
  paste0(s$prior, if (!is.na(s$lb)) paste0(" lb=", s$lb),
         if (!is.na(s$ub)) paste0(" ub=", s$ub), " [",
         paste(names(f), f, sep = "=", collapse = ", "), "]")
}
mk <- function(s) {
  a <- list(s$prior, class = s$class, coef = s$coef, group = s$group,
            resp = s$resp, dpar = s$dpar, nlpar = s$nlpar)
  if (!is.na(s$lb)) a$lb <- s$lb
  if (!is.na(s$ub)) a$ub <- s$ub
  do.call(ns$set_prior, a)
}
rowkey <- function(df) {
  df <- as.data.frame(df)
  for (c in c("class", "coef", "group", "resp", "dpar", "nlpar")) {
    if (is.null(df[[c]])) df[[c]] <- ""
    df[[c]][is.na(df[[c]])] <- ""
  }
  paste(df$class, df$coef, df$group, df$resp, df$dpar, df$nlpar, sep = "|")
}
reached <- function(design, pl) {
  r <- ns$resolve_priorlist(design, pl)
  pt <- design$frame[["par_template"]]
  nm <- vapply(r$entries, function(e) {
    nms <- ns$par_template_names(pt[[e$comp]], e$comp)
    paste(nms[e$idx], collapse = "+")
  }, "")
  list(txt = c(if (length(nm)) paste(nm, collapse = " ") else "(no density)",
               if (length(r$lower)) paste0("lb:", paste(names(r$lower), collapse = ",")),
               if (length(r$upper)) paste0("ub:", paste(names(r$upper), collapse = ","))),
       obj = r)
}
keep <- list()
for (case in names(r2_cases)) {
  cc <- r2_cases[[case]]
  m <- cc$build(ns_name)
  cat("\n== ", case, ": ", cc$what, " ==\n", sep = "")
  dp <- tryCatch(as.data.frame(ns$default_prior(m[[1]], data = d, family = m[[2]])),
                 error = function(e) {cat("MODEL REFUSED:", flat(conditionMessage(e)), "\n"); NULL})
  if (is.null(dp)) next
  dp <- dp[!dp$class %in% c("theta"), ]
  cat("-- default_prior rows (class|coef|group|resp|dpar|nlpar) --\n")
  cat(paste0("  ", sort(rowkey(dp))), sep = "\n")
  design <- if (who != "brms") ns$prior_design(m[[1]], d, m[[2]], list())
  for (i in seq_along(cc$specs)) {
    s <- cc$specs[[i]]
    r <- tryCatch({
      pl <- mk(s)
      vp <- as.data.frame(ns$validate_prior(pl, m[[1]], data = d, family = m[[2]]))
      u <- sort(rowkey(vp[vp$source == "user", ]))
      rc <- if (who != "brms") {
        rr <- reached(design, pl)
        keep[[paste(case, i)]] <- rr$obj
        paste(" || reach:", paste(rr$txt, collapse = " "))
      } else ""
      paste0("ACCEPT user rows: ", paste(u, collapse = " ; "), rc)
    }, error = function(e) paste0("REFUSE: ", substr(flat(conditionMessage(e)), 1, 400)))
    cat(sprintf("%-62s %s\n", lab(s), r))
  }
}
if (who != "brms") saveRDS(keep, file.path(R2_ROOT, "dev/mvprior-review2",
                                           paste0("r2-probe-", who, ".rds")))
