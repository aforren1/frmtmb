# Compare frmtmb's reading of every `||` case with brms 2.23.0's
# (lane dblbar): Z column for column, sd and cor names, the sd and cor
# rows of default_prior(), and whether the block correlates.
# Reads dev/dblbar-log/brms-probe.rds and frm-probe-<arm>.rds.
# Usage: Rscript dev/dblbar-compare.R base|lane
# Output: dev/dblbar-log/compare-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
dd <- dblbar_data()
br <- readRDS(file.path(wt, "dev/dblbar-log/brms-probe.rds"))
fr <- readRDS(file.path(wt, sprintf("dev/dblbar-log/frm-probe-%s.rds", arm)))
# brms's own renaming of a design column, as its coef names use it
rename <- function(x) {
  brms:::rename(sub("^[(]Intercept[)]$", "Intercept", x))
}
gkey <- function(g) {
  if (grepl("^mm[(]", g)) {
    return(paste0("mm", gsub("[^A-Za-z0-9_.,]", "", sub("^mm", "", g)) |>
                    gsub(pattern = ",", replacement = "")))
  }
  paste(sort(strsplit(g, ":", fixed = TRUE)[[1L]]), collapse = ":")
}
# brms's prefix of a row: resp, dpar (not mu) and nlpar joined by "_"
bpfx <- function(resp, dpar, nlpar) {
  p <- c(resp, if (!dpar %in% c("", "mu")) dpar, nlpar)
  paste(p[nzchar(p)], collapse = "_")
}
fpfx <- function(lp_key, multi, nlpars) {
  rd <- strsplit(lp_key, ".", fixed = TRUE)[[1L]]
  resp <- if (multi) rd[1L] else ""
  if (rd[2L] %in% nlpars) bpfx(resp, "", rd[2L]) else bpfx(resp, rd[2L], "")
}
brms_mats <- function(b, n) {
  rf <- b$frame_re
  out <- list()
  for (i in seq_len(nrow(rf))) {
    if (nzchar(rf$type[i])) next
    p <- bpfx(rf$resp[i], rf$dpar[i], rf$nlpar[i])
    zn <- paste0("Z_", rf$id[i], if (nzchar(p)) paste0("_", p), "_",
                 rf$cn[i])
    key <- paste(gkey(rf$group[i]), p, sep = "|")
    if (!is.null(b$Z[[zn]])) {
      # a multivariate model indexes its J by response
      J <- b$J[[paste0("J_", rf$id[i])]] %||%
        b$J[[paste0("J_", rf$id[i], "_", rf$resp[i])]]
      M <- matrix(0, n, max(J))
      M[cbind(seq_len(n), J)] <- b$Z[[zn]]
    } else {
      m <- 1L
      M <- NULL
      while (!is.null(b$Z[[paste0(zn, "_", m)]])) {
        J <- b$J[[paste0("J_", rf$id[i], "_", m)]]
        W <- b$W[[paste0("W_", rf$id[i], "_", m)]]
        if (is.null(M)) {
          L <- max(unlist(b$J[grep(paste0("^J_", rf$id[i], "_"),
                                   names(b$J))]))
          M <- matrix(0, n, L)
        }
        add <- matrix(0, n, ncol(M))
        add[cbind(seq_len(n), J)] <- W * b$Z[[paste0(zn, "_", m)]]
        M <- M + add
        m <- m + 1L
      }
      if (is.null(M)) next
    }
    out[[key]][[rf$coef[i]]] <- M
  }
  out
}
frm_mats <- function(f, multi, nlpars) {
  out <- list()
  for (bk in f$blocks) {
    for (cp in bk$comps) {
      key <- paste(gkey(bk$group), fpfx(cp$lp_key, multi, nlpars), sep = "|")
      for (j in seq_along(cp$cnms)) {
        nm <- rename(cp$cnms[j])
        old <- out[[key]][[nm]]
        out[[key]][[nm]] <- if (is.null(old)) cp$mats[[j]] else
          cbind(old, cp$mats[[j]])
      }
    }
  }
  out
}
# identical up to the order of the levels: every frmtmb level column,
# stacked over the coefficients, equals exactly one brms level column
same_up_to_levels <- function(Fl, Bl) {
  cf <- sort(names(Fl))
  if (!identical(cf, sort(names(Bl)))) return(FALSE)
  Fs <- do.call(rbind, Fl[cf])
  Bs <- do.call(rbind, Bl[cf])
  Bs <- Bs[, colSums(abs(Bs)) > 0, drop = FALSE]
  if (ncol(Fs) != ncol(Bs)) return(FALSE)
  m <- vapply(seq_len(ncol(Fs)), function(j) {
    hit <- which(colSums(abs(Bs - Fs[, j])) <= 1e-12 * max(1, abs(Fs[, j])))
    if (length(hit) == 1L) hit else NA_integer_
  }, 1L)
  !anyNA(m) && !anyDuplicated(m)
}
sink(file.path(wt, sprintf("dev/dblbar-log/compare-%s.txt", arm)),
     split = TRUE)
rows <- list()
for (nm in names(br)) {
  b <- br[[nm]]
  f <- fr[[nm]]
  berr <- b$error %||% ""
  ferr <- f$error %||% f$fit_error %||% ""
  if (nzchar(berr) || nzchar(ferr)) {
    rows[[nm]] <- data.frame(case = nm, brms_sds = NA, frm_sds = NA,
                             Z = NA, names = NA, prior = NA, cor = NA,
                             note = paste0(if (nzchar(berr)) "brms: error",
                                           if (nzchar(berr) && nzchar(ferr))
                                             "; ",
                                           if (nzchar(ferr)) "frmtmb: error"))
    cat(sprintf("\n%-7s brms: %s\n        frm:  %s\n", nm,
                substr(berr, 1, 90), substr(ferr, 1, 90)))
    next
  }
  multi <- grepl("^mvbf", b$case[1])
  nlpars <- unique(b$frame_re$nlpar[nzchar(b$frame_re$nlpar)])
  rf <- b$frame_re
  n <- length(b$Z[[1]])
  if (is.null(n) || !n) n <- 360L
  Bm <- brms_mats(b, n)
  Fm <- frm_mats(f, multi, nlpars)
  # every block of either side has its twin, with the same coefficients;
  # brms's sp (mo) coefficients have no Z in standata and are compared by
  # name only, so they are left out of frmtmb's side here
  sp <- rename(rf$coef[nzchar(rf$type)])
  Fm <- lapply(Fm, function(l) l[setdiff(names(l), sp)])
  Fm <- Fm[lengths(Fm) > 0L]
  zok <- length(Bm) > 0 && setequal(names(Bm), names(Fm)) &&
    all(vapply(names(Bm), function(k) {
      same_up_to_levels(Fm[[k]], Bm[[k]])
    }, NA))
  # sp (mo) coefficients have no Z in standata; their columns are
  # compared by name only
  # brms's sd and cor names, rebuilt from its frame_re
  bn <- character(0)
  for (id in unique(rf$id)) {
    r <- rf[rf$id == id, ]
    p <- mapply(bpfx, r$resp, r$dpar, r$nlpar)
    tn <- paste0(ifelse(nzchar(p), paste0(p, "_"), ""), r$coef)
    by <- r$by[1] %||% ""
    if (nzchar(by)) {
      # brms's get_rnames(): by-level first, the coefficients within it
      bl <- levels(dd[[by]])
      tn <- paste0(rep(tn, length(bl)), ":", by, rep(bl, each = length(tn)))
    }
    g <- r$group[1]
    bn <- c(bn, paste0("sd_", g, "__", tn))
    if (any(r$cor) && nrow(r) > 1L) {
      pr <- utils::combn(tn, 2)
      bn <- c(bn, paste0("cor_", g, "__", pr[1, ], "__", pr[2, ]))
    }
  }
  names_ok <- setequal(bn, f$names)
  pk <- function(p) {
    p <- p[!nzchar(p$coef), ]
    sort(paste(p$class, p$group, p$resp, p$dpar, p$nlpar, sep = "|"))
  }
  prior_ok <- identical(pk(b$prior), pk(f$prior))
  bcor <- any(rf$cor & ave(rf$id, rf$id, FUN = length) > 1)
  fcor <- any(grepl("^cor_", f$names))
  n_b_sd <- sum(grepl("^sd_", bn))
  n_f_sd <- sum(grepl("^sd_", f$names))
  rows[[nm]] <- data.frame(case = nm, brms_sds = n_b_sd, frm_sds = n_f_sd,
                           Z = zok, names = names_ok, prior = prior_ok,
                           cor = identical(bcor, fcor), note = "")
  cat(sprintf("\n%-7s sds brms %d frm %d | Z %s | names %s | prior %s | cor %s\n",
              nm, n_b_sd, n_f_sd, zok, names_ok, prior_ok,
              identical(bcor, fcor)))
  if (!names_ok) {
    cat("  brms only: ", paste(setdiff(bn, f$names), collapse = " "), "\n")
    cat("  frm only:  ", paste(setdiff(f$names, bn), collapse = " "), "\n")
  }
  if (!prior_ok) {
    cat("  brms rows: ", paste(pk(b$prior), collapse = "  "), "\n")
    cat("  frm rows:  ", paste(pk(f$prior), collapse = "  "), "\n")
  }
}
tab <- do.call(rbind, rows)
cat("\n")
print(tab, row.names = FALSE)
sink()
saveRDS(tab, file.path(wt, sprintf("dev/dblbar-log/compare-%s.rds", arm)))
