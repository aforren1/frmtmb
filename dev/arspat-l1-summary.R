# Lane arspat: the L1 and L4 numbers of dev/arspat-out/l1b-all.txt and
# l1c-sar.txt as one markdown table, for dev/arspat-findings.md.
#   Rscript dev/arspat-l1-summary.R > dev/arspat-out/l1-summary.md
rd <- function(f) readLines(f, warn = FALSE)
parse_log <- function(x, const_sar = 0) {
  cases <- list()
  cur <- NULL
  for (ln in x) {
    if (startsWith(ln, "######## ")) {
      cur <- sub("^######## ", "", ln)
      cases[[cur]] <- list(l1 = numeric(), grad = NA, l4 = NA, kind = "")
      next
    }
    if (is.null(cur)) next
    if (grepl("^L1 const at estimate", ln)) {
      v <- as.numeric(regmatches(ln, gregexpr("-?[0-9.]+e[-+][0-9]+", ln))[[1]])
      cases[[cur]]$l1 <- c(cases[[cur]]$l1, v[1])
      cases[[cur]]$grad <- v[2]
      cases[[cur]]$kind <- "joint"
    } else if (grepl("^L1 const at point", ln)) {
      v <- as.numeric(sub(".* ", "", ln))
      cases[[cur]]$l1 <- c(cases[[cur]]$l1, v)
    } else if (grepl("^marginal point", ln)) {
      # the printed const keeps every digit; where a constant is known
      # (brms's SAR densities leave out 0.5 N log(2 pi)) it is taken off
      # the two values, which are printed to 1e-10
      v <- if (const_sar == 0) {
        as.numeric(sub(".* const ", "", ln))
      } else {
        as.numeric(sub(".*stan ([-0-9.]+) .*", "\\1", ln)) -
          as.numeric(sub(".*ours ([-0-9.]+) .*", "\\1", ln)) - const_sar
      }
      cases[[cur]]$l1 <- c(cases[[cur]]$l1, v)
      cases[[cur]]$kind <- "marginal over Ymi"
    } else if (grepl("^L4 ", ln)) {
      cases[[cur]]$l4 <- as.numeric(sub(".* diff ([-0-9.e+]+) .*", "\\1", ln))
    } else if (grepl("^ERROR", ln)) {
      cases[[cur]]$kind <- paste("ERROR:", sub("^ERROR: *", "", ln))
    }
  }
  cases
}
a <- parse_log(rd("dev/arspat-out/l1b-all.txt"))
a <- a[!grepl("^mi_sar", names(a))]
b <- parse_log(rd("dev/arspat-out/l1c-sar.txt"), const_sar = 0.5 * 72 * log(2 * pi))
# fix round (dev/arspat-l1d.R): student registered terms, whose marginal
# is brms's Laplace value less c(nu, n, m) at nu held at 5, 10 and 30,
# where brms's student SAR densities leave out 0.5 N log(pi); and
# skew_normal's latent residual beside se()
s1 <- parse_log(rd("dev/arspat-out/l1d-student.txt"))
s1 <- s1[!grepl("^mi_sar", names(s1))]
s2 <- parse_log(rd("dev/arspat-out/l1d-student.txt"),
                const_sar = 0.5 * 72 * log(pi))
s2 <- s2[grepl("^mi_sar", names(s2))]
k <- parse_log(rd("dev/arspat-out/l1d-fix.txt"))
k <- k[grepl("^skew", names(k))]
all <- c(a, b, s1, s2, k)
cat("| case | identity | max abs const over 3 points | brms gradient at estimate | L4 brms max - frmtmb max |\n")
cat("|---|---|---|---|---|\n")
for (nm in names(all)) {
  r <- all[[nm]]
  cat(sprintf("| %s | %s | %s | %s | %s |\n", nm, r$kind,
              if (length(r$l1)) sprintf("%.1e (%d points)", max(abs(r$l1)),
                                        length(r$l1)) else "-",
              if (is.na(r$grad)) "-" else sprintf("%.1e", r$grad),
              if (is.na(r$l4)) "-" else sprintf("%.1e", r$l4)))
}
j <- Filter(function(r) r$kind == "joint", all)
m <- Filter(function(r) r$kind == "marginal over Ymi", all)
cat(sprintf("\nJoint cases: %d; max |const| %.1e; max |L4| %.1e over %d.\n",
            length(j), max(vapply(j, function(r) max(abs(r$l1)), 0)),
            max(abs(unlist(lapply(j, `[[`, "l4"))), na.rm = TRUE),
            sum(!is.na(unlist(lapply(j, `[[`, "l4"))))))
cat(sprintf("Marginal cases: %d; max |const| %.1e.\n", length(m),
            max(vapply(m, function(r) max(abs(r$l1)), 0))))
