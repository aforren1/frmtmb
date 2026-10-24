# Release review of 0.74.0, item 3: what a brms user meets, brms-style
# code with cox, car, spde and wiener, in a fresh process.
#   Rscript dev/review-rel074-user.R <libs;..> <attach;..|none> [brms]
# With "brms", brms is attached AFTER frmtmb (so brms::cox(), brms::bf()
# and brms::car() mask frmtmb's), as a brms user's script would.
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]], user))
if (isTRUE(a[3] == "brmsfirst")) suppressMessages(library(brms))
suppressMessages(library(frmtmb))
att <- strsplit(a[2], ";", fixed = TRUE)[[1]]
for (p in setdiff(att, "none")) suppressMessages(library(p, character.only = TRUE))
if (isTRUE(a[3] == "brms")) suppressMessages(library(brms))
inst <- function(p) nzchar(system.file(package = p))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "| installed:",
    paste(Filter(inst, c("frmtmb.survival", "frmtmb.spatial", "frmtmb.eam",
                         "frmtmb.latent")), collapse = " "),
    "| attached:", a[2], "| brms attached:", isTRUE(a[3] == "brms"), "\n")
g <- expand.grid(r = 1:4, c = 1:4)
n <- nrow(g)
W <- matrix(0, n, n)
for (i in 1:n) for (j in 1:n) {
  if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
}
dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
set.seed(20261023)
dd <- data.frame(loc = factor(rep(paste0("L", 1:n), each = 6),
                              levels = paste0("L", 1:n)),
                 x = rnorm(6 * n))
u <- rnorm(n, 0, 0.7)
dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$loc)] + rnorm(nrow(dd), 0, 0.5)
dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x))
dd$cen <- rbinom(nrow(dd), 1, 0.2)
nn <- 16; h <- 0.5
C0 <- diag(rep(h, nn)); C0[1, 1] <- C0[nn, nn] <- h / 2
G <- matrix(0, nn, nn)
for (i in seq_len(nn - 1)) {
  G[i, i] <- G[i, i] + 1 / h; G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
  G[i, i + 1] <- G[i + 1, i] <- -1 / h
}
fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
dd$node <- factor(as.integer(dd$loc))
set.seed(3)
dw <- data.frame(x = rnorm(120))
dw$rt <- 0.3 + rgamma(120, 2, 6)
dw$resp <- rbinom(120, 1, 0.6)
probe <- function(label, expr) {
  r <- tryCatch({
    f <- withCallingHandlers(eval(expr, envir = globalenv()),
                             warning = function(w) invokeRestart("muffleWarning"))
    if (inherits(f, "frmtmb_fit")) {
      paste("FIT logLik", format(as.numeric(logLik(f)), digits = 12))
    } else paste("VALUE", class(f)[1])
  }, error = function(e) {
    k <- class(e)[1]
    paste0("ERROR <", k, "> ", gsub("\\s+", " ", conditionMessage(e)))
  })
  cat(sprintf("%-34s %s\n", label, r))
}
assign("dd", dd, globalenv()); assign("W", W, globalenv())
assign("fem", fem, globalenv()); assign("dw", dw, globalenv())
probe("cox: family = cox()",
      quote(frm(bf(time | cens(cen) ~ x), family = cox(), data = dd)))
probe("cox: family = \"cox\"",
      quote(frm(bf(time | cens(cen) ~ x), family = "cox", data = dd)))
probe("cox: brmsfamily(\"cox\")",
      quote(frm(bf(time | cens(cen) ~ x), family = brmsfamily("cox"),
                data = dd)))
probe("cox: bf() + cox()",
      quote(frm(bf(time | cens(cen) ~ x) + cox(), data = dd)))
probe("car: car(W, gr = loc)",
      quote(frm(bf(y ~ x + car(W, gr = loc)), family = gaussian(),
                data = dd, data2 = list(W = W))))
probe("car: car(W, gr = loc, type=icar)",
      quote(frm(bf(y ~ x + car(W, gr = loc, type = "icar")),
                family = gaussian(), data = dd, data2 = list(W = W))))
probe("car: car(W) no gr",
      quote(frm(bf(y ~ x + car(W)), family = gaussian(), data = dd,
                data2 = list(W = W))))
probe("spde: spde(fem, gr = node)",
      quote(frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(),
                data = dd)))
probe("sar: sar(W)",
      quote(frm(bf(y ~ x + sar(W)), family = gaussian(), data = dd,
                data2 = list(W = W))))
probe("wiener: family = wiener()",
      quote(frm(bf(rt | dec(resp) ~ x), family = wiener(), data = dw)))
probe("wiener: family = \"wiener\"",
      quote(frm(bf(rt | dec(resp) ~ x), family = "wiener", data = dw)))
probe("wiener: brmsfamily(\"wiener\")",
      quote(frm(bf(rt | dec(resp) ~ x), family = brmsfamily("wiener"),
                data = dw)))
probe("dec() with gaussian",
      quote(frm(bf(rt | dec(resp) ~ x), family = gaussian(), data = dw)))
probe("bhaz() with cox",
      quote(frm(bf(time | cens(cen) + bhaz(df = 5) ~ x), family = "cox",
                data = dd)))
probe("mixture_mvn",
      quote(frm(bf(cbind(y, x) ~ 1), family = mixture_mvn(K = 2, D = 2),
                data = dd)))
probe("royston_parmar",
      quote(frm(bf(time | cens(cen) ~ x), family = royston_parmar(),
                data = dd)))
