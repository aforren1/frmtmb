# The two pairs against the baseline. Log dev/impsamp-out/equiv-verdicts.txt
b <- read.delim("C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out/results.tsv")
a <- read.delim("dev/impsamp-equiv-out/results.tsv")
cols <- c("id", "verdict", "verdict_measured", "level", "frm_ok", "l1_const", "l1_grad", "l1_ok", "l4_ok", "l4_z", "l4_control", "seconds")
print(a[a$id %in% c("p124", "p069"), cols], row.names = FALSE)
m <- merge(b[, c("id", "verdict_measured")], a[, c("id", "verdict_measured")], by = "id", suffixes = c(".b", ".a"))
print(m[m$verdict_measured.b != m$verdict_measured.a, ], row.names = FALSE)
