# Counts from the off-grid run, generated. Log dev/impsamp-out/offgrid-summary.txt
root <- "dev/impsamp-out/offgrid"
rd <- function(d) {
  f <- list.files(file.path(root, d), "^[^-]*-[0-9]+[.]tsv$", full.names = TRUE)
  do.call(rbind, lapply(f, function(x) {
    if (grepl("\tERROR\t", readLines(x, 1L))) {
      cat("ERROR in", x, "\n")
      return(NULL)
    }
    read.delim(x, quote = "", colClasses = c(warn = "character"))
  }))
}
com <- c("design", "seed", "rounds", "n_warn", "max_z_est", "max_z_lap", "newton",
         "se_ratio_min", "se_ratio_max", "ess_min")
a1 <- rd("accuracy"); a2 <- rd("accuracy2"); m <- rd("multid")
scale_d <- c("scale1e5", "pois1e5", "year")
all <- rbind(a1[com], a2[!a2$design %in% scale_d, com], m[com])
cat("fits", nrow(all), "designs", length(unique(all$design)), "\n")
print(table(all$rounds))
bad <- all$design == "re_sigma_small" & all$seed == 1
cat("excl re_sigma_small 1: max z", max(all$max_z_est[!bad], na.rm = TRUE),
    "newton", max(all$newton[!bad], na.rm = TRUE),
    "lap max z", max(all$max_z_lap[!bad], na.rm = TRUE), "\n")
cat("SE ratio range", range(all$se_ratio_min[!bad], na.rm = TRUE),
    range(all$se_ratio_max[!bad], na.rm = TRUE), "\n")
print(all[which(all$se_ratio_max > 1.08 | all$se_ratio_min < 0.92),
          c("design", "seed", "se_ratio_min", "se_ratio_max", "n_warn", "ess_min")])
cat("warned", sum(all$n_warn > 0), "\n")
print(all[all$n_warn > 0 & all$ess_min >= 0.25, c("design", "seed", "n_warn", "ess_min")])
cat("the scale designs of accuracy2:\n")
print(a2[a2$design %in% scale_d, com], row.names = FALSE)
