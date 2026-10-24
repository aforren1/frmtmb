# frmtmb.sample on car() and spde() fits, before and after the move: the
# same seed must give the same draws, and the block must sample
# centered, with the reason frmtmb.sample names.
#   Rscript dev/spatial-sample-check.R <out.rds>   (SPATIAL_LIBS, SPATIAL_PKGS)
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
pk <- strsplit(Sys.getenv("SPATIAL_PKGS"), ",", fixed = TRUE)[[1]]
for (p in pk) suppressMessages(library(p, character.only = TRUE))
g <- expand.grid(r = 1:3, c = 3)
W <- matrix(0, 6, 6, dimnames = list(paste0("L", 1:6), paste0("L", 1:6)))
W[cbind(1:5, 2:6)] <- 1
W <- W + t(W)
set.seed(2)
d <- data.frame(loc = factor(rep(rownames(W), each = 8), levels = rownames(W)),
                x = rnorm(48))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$loc] + rnorm(48, 0, 0.5)
fc <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")), family = gaussian(),
          data = d, data2 = list(W = W))
nn <- 8; h <- 0.5
C0 <- diag(rep(h, nn)); C0[1, 1] <- h / 2; C0[nn, nn] <- h / 2
G <- matrix(0, nn, nn)
for (i in seq_len(nn - 1)) {
  G[i, i] <- G[i, i] + 1 / h; G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
  G[i, i + 1] <- -1 / h; G[i + 1, i] <- -1 / h
}
fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
ds <- data.frame(node = factor(rep(seq_len(nn), each = 6)), x = rnorm(48))
ds$y <- 0.3 + 0.4 * ds$x + rnorm(nn)[ds$node] + rnorm(48, 0, 0.3)
fs <- frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(), data = ds)
out <- list()
for (nm in c("car", "spde")) {
  fit <- if (nm == "car") fc else fs
  msgs <- character()
  dr <- withCallingHandlers(
    frm_sample(fit, chains = 1, iter = 300, warmup = 150, seed = 11,
               refresh = 0),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  out[[nm]] <- list(draws = as.matrix(as_draws_matrix(dr)), msgs = msgs)
  cat(nm, "draws", dim(out[[nm]]$draws), "\n")
  cat(paste(" msg:", substr(msgs, 1, 160)), sep = "\n")
}
saveRDS(out, commandArgs(trailingOnly = TRUE)[1])
