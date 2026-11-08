# Lane bigx: the survey's metrics, lane against the record.
rec <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
new <- do.call(rbind, lapply(list.files(Sys.getenv("BIGX_EQUIV_OUT", "dev/bigx-equiv-out"), "[.]tsv$", full.names = TRUE), function(f) read.delim(f, stringsAsFactors = FALSE)[, c("id", "l1_const", "l1_grad", "l2_metric", "l4_z", "frm_ok", "brms_ok")]))
m <- merge(rec[, c("id", "l1_const", "l1_grad", "l2_metric", "l4_z", "frm_ok", "brms_ok")], new, by = "id")
d <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a - b) / pmax(1e-8, abs(a)))
for (v in c("l1_const", "l1_grad", "l2_metric")) {
  x <- d(m[[paste0(v, ".x")]], m[[paste0(v, ".y")]])
  cat(v, ": pairs differing by more than 1e-6 relative:", sum(x > 1e-6, na.rm = TRUE), "; NA on one side:", sum(is.na(x)), "\n")
  if (any(x > 1e-6, na.rm = TRUE)) print(m[which(x > 1e-6), c("id", paste0(v, ".x"), paste0(v, ".y"))])
}
cat("l4_z text differs:", sum(m$l4_z.x %in% NA != m$l4_z.y %in% c(NA, "") | (!is.na(m$l4_z.x) & m$l4_z.x != m$l4_z.y)), "\n")
cat("frm_ok differs:", sum(m$frm_ok.x != m$frm_ok.y, na.rm = TRUE), " brms_ok differs:", sum(m$brms_ok.x != m$brms_ok.y, na.rm = TRUE), "\n")