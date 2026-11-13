# Reviewer of lane priortab: `+` duplicates and order, response names,
# extra tables, against brms 2.23.0 (no compile).
# Rscript dev/reviews/priortab-review-probe1.R > dev/reviews/priortab-review-probe1.txt
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  loadNamespace("brms"); loadNamespace("frmtmb"); loadNamespace("frmtmb.sample")
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 160, warn = 1)
B <- asNamespace("brms"); F <- asNamespace("frmtmb")
cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "lb",
          "ub", "source")
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
         }, message = function(m) {
           cat("MSG: ", conditionMessage(m)); invokeRestart("muffleMessage")
         }),
       error = function(e) {
         cat("ERROR: ", conditionMessage(e), "\n")
         NULL
       })
  if (is.data.frame(r)) {
    r <- as.data.frame(r)
    print(r[, intersect(cols, names(r)), drop = FALSE], row.names = TRUE)
  } else if (!is.null(r)) {
    print(r)
  }
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$Y_a <- d$y
d$Ya <- d$y2
d$y.1 <- d$y
d$y.2 <- d$y2

# --- `+` with a duplicated slot, both orders, class b and coef z
for (pk in c("brms", "frmtmb")) {
  P <- asNamespace(pk)
  tab <- P$default_prior(y ~ z + x, data = d)
  sp_b <- P$set_prior("normal(0, 1)", class = "b")
  sp_z <- P$set_prior("normal(0, 1)", class = "b", coef = "z")
  sp_s <- P$set_prior("normal(0, 1)", class = "sigma")
  for (nm in c("sp_b", "sp_z", "sp_s")) {
    sp <- get(nm)
    show(paste(pk, "tab +", nm, "validate"),
         P$validate_prior(tab + sp, y ~ z + x, data = d))
    show(paste(pk, nm, "+ tab validate"),
         P$validate_prior(sp + tab, y ~ z + x, data = d))
  }
  # an edited table + a set_prior on another slot: which wins, row order
  tab2 <- tab
  tab2$prior[tab2$class == "b" & tab2$coef == "z"] <- "normal(0, 5)"
  show(paste(pk, "edited tab + sp_s"),
       P$validate_prior(tab2 + sp_s, y ~ z + x, data = d))
  show(paste(pk, "tab + tab"), P$validate_prior(tab + tab, y ~ z + x,
                                                data = d))
}
# brms: what does + of two set_prior on the same slot do
show("brms sp + sp same slot", B$validate_prior(
  B$set_prior("normal(0, 1)", class = "b") +
    B$set_prior("normal(0, 5)", class = "b"), y ~ z, data = d))
show("frmtmb sp + sp same slot", F$validate_prior(
  F$set_prior("normal(0, 1)", class = "b") +
    F$set_prior("normal(0, 5)", class = "b"), y ~ z, data = d))
# frmtmb sampling-route table + sp: route kept
tabs <- F$default_prior(y ~ z + x, data = d, route = "sample")
show("frmtmb sample tab + sp coef z, route attr", {
  r <- tabs + F$set_prior("normal(0, 1)", class = "b", coef = "x")
  attr(r, "route")
})
show("frmtmb fit tab + sample tab", {
  r <- F$default_prior(y ~ z + x, data = d) + tabs
  attr(r, "route")
})

# --- response names
show("brms mvbf(Y_a, Ya)", B$default_prior(
  B$mvbf(B$bf(Y_a ~ z), B$bf(Ya ~ z), rescor = FALSE), data = d))
show("frmtmb mvbf(Y_a, Ya)", F$default_prior(
  F$mvbf(F$bf(Y_a ~ z), F$bf(Ya ~ z), rescor = FALSE), data = d))
show("frmtmb frm mvbf(Y_a, Ya) resp = Ya", {
  f <- F$frm(F$mvbf(F$bf(Y_a ~ z), F$bf(Ya ~ z), rescor = FALSE), data = d,
             prior = F$set_prior("normal(5, 0.01)", class = "b", resp = "Ya"))
  F$fixef(f)
})
show("brms mvbind(y.1, y.2)", B$default_prior(
  B$bf(mvbind(y.1, y.2) ~ z) + B$set_rescor(FALSE), data = d))
show("frmtmb mvbind(y.1, y.2)", F$default_prior(
  F$bf(mvbind(y.1, y.2) ~ z) + F$set_rescor(FALSE), data = d))
show("frmtmb validate resp = y.1", F$validate_prior(
  F$set_prior("normal(0, 1)", class = "b", resp = "y.1"),
  F$bf(mvbind(y.1, y.2) ~ z) + F$set_rescor(FALSE), data = d))
show("frmtmb validate resp = y1", F$validate_prior(
  F$set_prior("normal(0, 1)", class = "b", resp = "y1"),
  F$bf(mvbind(y.1, y.2) ~ z) + F$set_rescor(FALSE), data = d))
show("brms validate resp = y.1", B$validate_prior(
  B$set_prior("normal(0, 1)", class = "b", resp = "y.1"),
  B$bf(mvbind(y.1, y.2) ~ z) + B$set_rescor(FALSE), data = d))
show("frmtmb frm table passed back y.1", {
  f0 <- F$bf(mvbind(y.1, y.2) ~ z) + F$set_rescor(FALSE)
  tb <- F$default_prior(f0, data = d)
  tb$prior[tb$class == "b" & tb$coef == "z" & tb$resp == "y1"] <-
    "normal(5, 0.01)"
  fit <- F$frm(f0, data = d, prior = tb)
  list(F$fixef(fit), as.data.frame(fit$prior))
})
show("frmtmb Y_a table + set_prior resp=Y_a (dup?)", {
  f0 <- F$bf(mvbind(Y_a, y2) ~ z) + F$set_rescor(FALSE)
  tb <- F$default_prior(f0, data = d)
  tb <- tb[tb$coef != "z", ]
  F$validate_prior(tb + F$set_prior("normal(0, 1)", class = "b",
                                    coef = "z", resp = "Y_a"), f0, data = d)
})
# --- cor with resp still accepted (Q2)
show("frmtmb cor resp accepted", F$validate_prior(
  F$set_prior("lkj(3)", class = "cor", resp = "y2"),
  F$bf(mvbind(y, y2) ~ z + (1 + z | g)) + F$set_rescor(FALSE), data = d))
# --- extra tables
show("brms brmsfamily gaussian identity sigma", B$default_prior(
  y ~ z, data = d, family = B$brmsfamily("gaussian", link_sigma = "identity")))
show("frmtmb brmsfamily gaussian identity sigma (fit)", F$default_prior(
  y ~ z, data = d, family = F$brmsfamily("gaussian", link_sigma = "identity")))
show("frmtmb brmsfamily gaussian identity sigma (sample)", F$default_prior(
  y ~ z, data = d, family = F$brmsfamily("gaussian", link_sigma = "identity"),
  route = "sample"))
show("brms student identity nu", B$default_prior(
  y ~ z, data = d, family = B$student(link_nu = "identity")))
show("frmtmb student identity nu fit", F$default_prior(
  y ~ z, data = d, family = F$student(link_nu = "identity")))
show("brms Beta identity phi", B$default_prior(
  pb ~ z, data = transform(d, pb = plogis(y)),
  family = B$Beta(link_phi = "identity")))
show("frmtmb Beta identity phi fit", F$default_prior(
  pb ~ z, data = transform(d, pb = plogis(y)),
  family = F$Beta(link_phi = "identity")))
show("brms nl lb in prior", B$validate_prior(
  B$set_prior("normal(1, 1)", nlpar = "a", lb = 0) +
    B$set_prior("normal(0, 1)", nlpar = "b"),
  B$bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, nl = TRUE), data = d))
show("frmtmb nl lb in prior", F$validate_prior(
  F$set_prior("normal(1, 1)", nlpar = "a", lb = 0) +
    F$set_prior("normal(0, 1)", nlpar = "b"),
  F$bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, nl = TRUE), data = d))
cat("\ndone\n")
