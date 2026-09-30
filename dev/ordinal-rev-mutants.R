# Reviewer, lane ordinal: build mutants of the lane's core source (a copy
# in the session scratchpad; the worktree is not touched) and install
# each into its own scratch library.
# Output: dev/ordinal-rev-log-mutants.txt
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
src0 <- file.path(sp, "ordrev-src", "core")
muts <- list(
  M1_delta_offset = list("R/thres.R",
    "return(r[1L] + (seq_len(k) - 1) * d)", "return(r[1L] + seq_len(k) * d)"),
  M1b_delta_sign = list("R/thres.R",
    "return(r[1L] + (seq_len(k) - 1) * d)",
    "return(r[1L] - (seq_len(k) - 1) * d)"),
  M2_acat_logit_nodisc = list("R/families.R",
    "      E <- E * disc\n", "      E <- E\n"),
  M2b_grouped_seq_nodisc = list("R/thres.R",
    "    Mj <- disc * (tau[s + pmin(j, nk) - 1L] - eta)",
    "    Mj <- (tau[s + pmin(j, nk) - 1L] - eta)"),
  M2c_cratio_nodisc = list("R/families.R",
    "      M <- M * disc\n", "      M <- M * (if (stopping) disc else 1)\n"),
  M3_stz_ordered_uncentered = list("R/thres.R",
    "return(u - sum(u) / k)", "return(u - u[1L])"),
  M3b_stz_inverse = list("R/thres.R",
    "sum_to_zero = if (ordered) log(diff(tg)) else tg[-k])",
    "sum_to_zero = if (ordered) log(diff(tg)) else tg[-1L])"),
  M3c_stz_half = list("R/thres.R", "  c(r, -sum(r))\n", "  c(r, -sum(r) / 2)\n"),
  M4_delta_value = list("R/confint.R",
    "value = if (ordered) exp else identity", "value = identity"),
  M5_delta_prior_scale = list("R/priors.R",
    "scale = if (ordered) \"sd\" else \"internal\", link = NULL,",
    "scale = \"internal\", link = NULL,"),
  M6_catprobs_nodisc = list("R/families.R", "  M <- disc * M\n", "  M <- M\n")
)
R <- file.path(R.home("bin"), "R.exe")
sink(file.path(wt, "dev/ordinal-rev-log-mutants.txt"), split = TRUE)
for (id in names(muts)) {
  m <- muts[[id]]
  dst <- file.path(sp, "ordrev-mut", id)
  lib <- file.path(sp, "ordrev-mutlib", id)
  unlink(dst, recursive = TRUE); unlink(lib, recursive = TRUE)
  dir.create(dst, recursive = TRUE); dir.create(lib, recursive = TRUE)
  file.copy(list.files(src0, full.names = TRUE), dst, recursive = TRUE)
  f <- file.path(dst, m[[1]])
  txt <- paste(readLines(f), collapse = "\n")
  hits <- lengths(regmatches(txt, gregexpr(m[[2]], txt, fixed = TRUE)))
  txt <- sub(m[[2]], m[[3]], txt, fixed = TRUE)
  writeLines(txt, f)
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", lib),
                     "--no-multiarch", "--no-test-load", dst),
                stdout = TRUE, stderr = TRUE,
                env = paste0("R_LIBS=", paste(c(lib,
                  "C:/Users/adf44/source/r/wt-ordinal-lib",
                  "C:/Users/adf44/source/r/rellib-r4",
                  "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
                  collapse = ";")))
  cat(id, ": matches", hits, " install:", tail(st, 1), "\n")
}
sink()
