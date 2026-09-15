# Summarize the reviewer's probe. Every column is TEST against the
# CONTROL measured in a separate process with the same load order and
# frmtmb.sample never loaded.
root <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out/probe"
modes <- c("S", "T", "U", "N", "P", "Q", "G", "D", "R",
           "RN", "RN2", "BR", "CA", "CA2")
label <- c(S = "brms, then frmtmb.sample",
           T = "frmtmb.sample, then brms",
           U = "frmtmb.sample; brms only loaded",
           N = "frmtmb.sample alone",
           P = "7 other owners, then frmtmb.sample",
           Q = "frmtmb.sample, then 7 other owners",
           G = "gratia, then frmtmb.sample",
           D = "brms, sample; detach and reattach both",
           R = "loo, sample; loo unloaded; brms loaded",
           RN = "brms attached; sample requireNamespace ONLY",
           RN2 = "sample requireNamespace ONLY; then brms",
           BR = "both attached; brms unloaded and RELOADED",
           CA = "brms, frmtmb, then frmtmb.sample",
           CA2 = "brms, frmtmb.sample, then frmtmb")

key <- function(r) paste0(r$via, "|", r$meth, "|", r$from)

cat("== reviewer collision probe, dev/sgrev-probe.R ==\n")
cat("one process per cell; TEST vs a CONTROL process with the same\n")
cat("order and frmtmb.sample absent. lookup = caller env, then the\n")
cat(".__S3MethodsTable__. of environment(generic).\n")
cat("brmsfit: control reaches a brmsfit method, test does not reach",
    "the same one\n")
cat("draws:   frmtmb_draws method that would run is not",
    "frmtmb.sample's (a foreign .default counts as lost)\n")
cat("owners:  (name,class) pairs the control resolves to a CLASS",
    "method and the test does not reach\n")
cat("gS3:     the same comparison through getS3method()\n\n")
cat(sprintf("%-5s %-4s %-42s %-9s %-7s %-9s %-9s\n",
            "build", "mode", "order", "brmsfit", "draws", "owners",
            "gS3"))

legs <- character()
for (build in c("BASE", "FIX")) {
  for (mode in modes) {
    ft <- file.path(root, paste0(build, "-", mode, "-test.rds"))
    fc <- file.path(root, paste0(build, "-", mode, "-control.rds"))
    if (!file.exists(ft) || !file.exists(fc)) next
    tt <- readRDS(ft); cc <- readRDS(fc)
    nms <- tt$base_defined
    bf_d <- 0; bf_n <- 0
    ow_d <- 0; ow_n <- 0
    g3_d <- 0; g3_n <- 0
    dr_d <- 0
    lost_names <- character(); lost_pairs <- character()
    for (nm in nms) {
      T1 <- tt$res[[nm]]; C1 <- cc$res[[nm]]
      # brmsfit
      cr <- C1$rows[["brmsfit"]]
      if (!is.null(cr) && identical(cr$meth, "brmsfit")) {
        bf_n <- bf_n + 1
        trw <- T1$rows[["brmsfit"]]
        if (is.null(trw) || !identical(key(trw), key(cr))) {
          bf_d <- bf_d + 1
          lost_names <- c(lost_names, nm)
        }
      }
      # owner class methods
      for (cl in names(C1$rows)) {
        if (cl %in% c("brmsfit", "frmtmb_draws")) next
        crw <- C1$rows[[cl]]
        if (is.null(crw) || !identical(crw$meth, cl)) next
        ow_n <- ow_n + 1
        trw <- T1$rows[[cl]]
        if (is.null(trw) || !identical(key(trw), key(crw))) {
          ow_d <- ow_d + 1
          lost_pairs <- c(lost_pairs, paste0(nm, ".", cl))
        }
      }
      # getS3method
      for (cl in names(C1$gs3)) {
        cv <- C1$gs3[[cl]]
        if (is.na(cv)) next
        if (cl %in% "frmtmb_draws") next
        g3_n <- g3_n + 1
        tv <- T1$gs3[[cl]]
        if (!identical(tv, cv)) g3_d <- g3_d + 1
      }
      # draws
      trw <- T1$rows[["frmtmb_draws"]]
      ours <- !is.null(trw) && identical(trw$meth, "frmtmb_draws") &&
        identical(trw$from, "frmtmb.sample")
      if (!ours) {
        dr_d <- dr_d + 1
        lost_pairs <- c(lost_pairs,
                        paste0("DRAWS:", nm, "<-",
                               if (is.null(trw)) "NA" else
                                 paste0(trw$from, ".", trw$meth)))
      }
      for (r in T1$rows) if (identical(r$via, "callerenv"))
        legs <- c(legs, paste(build, mode, nm, r$class))
    }
    cat(sprintf("%-5s %-4s %-42s %-9s %-7s %-9s %-9s\n",
                build, mode, label[[mode]],
                paste0(bf_d, "/", bf_n),
                paste0(dr_d, "/", length(nms)),
                paste0(ow_d, "/", ow_n),
                paste0(g3_d, "/", g3_n)))
    if (length(lost_pairs)) {
      lp <- grep("^DRAWS:", lost_pairs, value = TRUE)
      if (length(lp))
        cat("      draws lost: ", paste(lp, collapse = " "), "\n")
    }
    if (length(tt$notes))
      cat("      notes: ", paste(tt$notes, collapse = " | "), "\n")
    act <- vapply(tt$res, function(r) isTRUE(r$active_ns), NA)
    cat("      active in namespace: ", sum(act, na.rm = TRUE), "/",
        length(act), "\n", sep = "")
  }
}
cat("\nresolutions that came from the CALLER-ENV leg (leg 1): ",
    length(legs), "\n")
if (length(legs)) cat(paste(head(unique(legs), 20), collapse = "\n"), "\n")
cat("DONE\n")
