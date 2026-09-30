# Reviewer, punch round 2: try to break cmc_changes_columns(). A plain bar
# is compared with brms's standata (Z_k Z_k', which does not depend on
# level order); a level-read structure is checked for its refusal.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
set.seed(61)
n <- 240
d <- data.frame(g = factor(rep(sprintf("G%02d", 1:24), each = 10)),
                f = factor(sample(c("a", "b", "c"), n, TRUE),
                           levels = c("a", "b", "c")),
                f2 = factor(sample(c("p", "q"), n, TRUE)),
                x = rnorm(n), xi = sample(1:3, n, TRUE))
d$fu <- factor(as.character(d$f), levels = c("a", "b", "c", "zz"))
d$fchr <- as.character(d$f)
d$o <- factor(as.character(d$f), levels = c("a", "b", "c"), ordered = TRUE)
d$fna <- d$f; d$fna[sample(n, 15)] <- NA
d$y <- rnorm(n, as.numeric(d$f) + d$x) + rnorm(24, 0, 0.5)[d$g]
cmpZ <- function(lhs, cmc = FALSE) {
  f <- as.formula(sprintf("y ~ x + (%s | g)", lhs))
  fr <- q(frm(bf(f, cmc = cmc), data = d, dry_run = "frame"))
  sd <- q(brms::standata(brms::bf(f, cmc = cmc), d))
  if (inherits(fr, "ERR") || inherits(sd, "ERR")) {
    return(sprintf("frm: %s | brms: %s",
                   if (inherits(fr, "ERR")) substr(fr, 1, 90) else "ok",
                   if (inherits(sd, "ERR")) substr(sd, 1, 90) else "ok"))
  }
  cn <- unlist(lapply(fr$re_blocks, function(b)
    lapply(b$components, `[[`, "cnms")))
  K <- length(cn); nn <- sd$N
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  kb <- length(grep("^Z_1_[0-9]+$", names(sd)))
  ok <- K == kb && nrow(Z) == nn
  if (ok) for (k in seq_len(K)) {
    Zf <- Z[, seq(k, ncol(Z), by = K), drop = FALSE]
    Zb <- matrix(0, nn, sd$N_1)
    Zb[cbind(seq_len(nn), sd$J_1)] <- sd[[paste0("Z_1_", k)]]
    ok <- ok && isTRUE(all.equal(unname(tcrossprod(Zf)),
                                 unname(tcrossprod(Zb))))
  }
  sprintf("n=%d cnms [%s] brms coefs %d, Z Z' equal %s", nrow(Z),
          paste(cn, collapse = ","), kb, ok)
}
cs_dec <- function(lhs) {
  f <- as.formula(sprintf("y ~ x + cs(%s | g)", lhs))
  a <- q(frm(bf(f), data = d, dry_run = "frame"))
  b <- q(frm(bf(f, cmc = FALSE), data = d, dry_run = "frame"))
  cn <- function(fr) if (inherits(fr, "ERR"))
    paste("REFUSED:", substr(fr, 1, 70)) else
      paste(unlist(lapply(fr$re_blocks, function(bk)
        lapply(bk$components, `[[`, "cnms"))), collapse = ",")
  sprintf("cs default [%s] | cs cmc=FALSE [%s]", cn(a), cn(b))
}
shapes <- c("0 + f:f2", "0 + f + f2", "0 + f * x", "0 + fu", "0 + fchr",
            "0 + o", "0 + factor(xi)", "0 + fna", "0 + f:x", "0 + x",
            "0 + fu:x", "0 + o:x")
for (s in shapes) {
  cat(sprintf("%-16s plain cmc=FALSE: %s\n", s, cmpZ(s, FALSE)))
  cat(sprintf("%-16s plain cmc=TRUE : %s\n", s, cmpZ(s, TRUE)))
  cat(sprintf("%-16s %s\n", s, cs_dec(s)))
}
cat("== direct calls of the helper on the model frame\n")
fr <- frm(bf(y ~ x + (0 + f | g) + (0 + fu | g) + (0 + o | g) +
               (0 + factor(xi) | g) + (0 + f:f2 | g) + (0 + fchr | g)),
          data = d, dry_run = "frame")
mf <- fr$data_frame %||% fr[["mf"]]
cat("frame has data_frame:", !is.null(fr$data_frame), "\n")
cat("ar1(0 + t | g) over a factor still refused:\n")
d$t <- factor(rep(1:10, 24))
cat(" ", substr(q(frm(bf(y ~ x + ar1(0 + t | g), cmc = FALSE), data = d,
                      dry_run = "frame")), 1, 90), "\n")
cat("DONE\n")
