# Side by side: brms verdict and user rows, lane, base. Also compare
# default_prior row sets per case, and the resolved entries lane vs base
# with identical().
root <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2"
rd <- function(w) readLines(file.path(root, paste0("r2-probe-", w, ".txt")))
parse_probe <- function(L) {
  out <- list(); rows <- list(); case <- NA; inrows <- FALSE
  for (l in L) {
    if (grepl("^== ", l)) { case <- sub("^== ([^:]+):.*", "\\1", l); inrows <- FALSE; next }
    if (grepl("^-- default_prior", l)) { inrows <- TRUE; rows[[case]] <- character(0); next }
    if (grepl("^MODEL REFUSED", l)) { rows[[case]] <- "MODEL REFUSED"; next }
    if (inrows && grepl("^  ", l)) { rows[[case]] <- c(rows[[case]], trimws(l)); next }
    if (!is.na(case) && grepl("(ACCEPT|REFUSE)", l)) {
      inrows <- FALSE
      k <- regexpr("(ACCEPT|REFUSE)", l)
      spec <- trimws(substr(l, 1, k - 1)); v <- substr(l, k, nchar(l))
      out[[length(out) + 1]] <- data.frame(case = case, spec = spec,
        verdict = substr(v, 1, 6),
        user = if (startsWith(v, "ACCEPT")) sub(" [|][|] reach:.*", "", sub("ACCEPT user rows: ", "", v)) else "",
        reach = if (grepl("reach:", v)) sub(".* [|][|] reach: ", "", v) else "",
        msg = if (startsWith(v, "REFUSE")) substr(v, 9, 200) else "")
    }
  }
  list(specs = do.call(rbind, out), rows = rows)
}
b <- parse_probe(rd("brms")); l <- parse_probe(rd("lane")); s <- parse_probe(rd("base"))
m <- merge(merge(b$specs, l$specs, by = c("case", "spec"), suffixes = c(".brms", ".lane"), all = TRUE),
           setNames(s$specs, c("case", "spec", paste0(names(s$specs)[-(1:2)], ".base"))),
           by = c("case", "spec"), all = TRUE)
m$agree_lane <- m$verdict.brms == m$verdict.lane
m$agree_base <- m$verdict.brms == m$verdict.base
cat("specs:", nrow(m), " lane agrees with brms:", sum(m$agree_lane, na.rm = TRUE),
    " base agrees:", sum(m$agree_base, na.rm = TRUE), "\n")
# user rows agreement where both accept
both <- m$verdict.brms == "ACCEPT" & m$verdict.lane == "ACCEPT"
ur <- function(x) gsub("[(]Intercept[)]", "Intercept", x)
m$rows_same <- ifelse(both, ur(m$user.brms) == ur(m$user.lane), NA)
cat("both accept:", sum(both, na.rm = TRUE), " user rows same:", sum(m$rows_same, na.rm = TRUE), "\n")
cat("\n--- verdict disagreements lane vs brms ---\n")
dd <- m[!m$agree_lane | is.na(m$agree_lane), ]
for (i in seq_len(nrow(dd))) cat(sprintf("[%s] %s\n   brms=%s lane=%s base=%s\n   lane msg: %s\n",
  dd$case[i], dd$spec[i], dd$verdict.brms[i], dd$verdict.lane[i], dd$verdict.base[i], dd$msg.lane[i]))
cat("\n--- both accept but user rows differ ---\n")
dd <- m[which(both & !m$rows_same), ]
for (i in seq_len(nrow(dd))) cat(sprintf("[%s] %s\n   brms: %s\n   lane: %s\n", dd$case[i], dd$spec[i], dd$user.brms[i], dd$user.lane[i]))
cat("\n--- accepted on lane: reach, and base reach ---\n")
dd <- m[which(m$verdict.lane == "ACCEPT"), ]
for (i in seq_len(nrow(dd))) cat(sprintf("[%s] %-50s lane: %s | base: %s%s\n", dd$case[i], dd$spec[i],
  dd$reach.lane[i], if (identical(dd$verdict.base[i], "ACCEPT")) dd$reach.base[i] else paste("REFUSED", substr(dd$msg.base[i], 1, 60)),
  if (identical(dd$reach.lane[i], dd$reach.base[i])) "  [same]" else "  [DIFF]"))
cat("\n--- verdict changes base -> lane ---\n")
dd <- m[which(m$verdict.lane != m$verdict.base), ]
for (i in seq_len(nrow(dd))) cat(sprintf("[%s] %s base=%s lane=%s brms=%s | base reach: %s\n", dd$case[i], dd$spec[i],
  dd$verdict.base[i], dd$verdict.lane[i], dd$verdict.brms[i], dd$reach.base[i]))
cat("\n--- default_prior rows, per case (brms vs lane vs base), class theta dropped ---\n")
norm <- function(r) sort(unique(gsub("[(]Intercept[)]", "Intercept", r)))
for (cs in names(b$rows)) {
  br <- norm(b$rows[[cs]]); lr <- norm(l$rows[[cs]]); sr <- norm(s$rows[[cs]])
  cat(sprintf("%-16s brms %2d lane %2d base %2d | lane-only: %s | brms-only: %s | lane==base: %s\n", cs,
    length(br), length(lr), length(sr), paste(setdiff(lr, br), collapse = " "),
    paste(setdiff(br, lr), collapse = " "), identical(lr, sr)))
}
# resolved entries identical lane vs base where both accept
kl <- readRDS(file.path(root, "r2-probe-lane.rds")); kb <- readRDS(file.path(root, "r2-probe-base.rds"))
com <- intersect(names(kl), names(kb))
same <- vapply(com, function(k) identical(kl[[k]], kb[[k]]), TRUE)
cat("\nresolved objects both arms accept:", length(com), " identical():", sum(same), "\n")
if (any(!same)) print(com[!same])
cat("lane-only accepted:", paste(setdiff(names(kl), names(kb)), collapse = ", "), "\n")
cat("base-only accepted:", paste(setdiff(names(kb), names(kl)), collapse = ", "), "\n")
