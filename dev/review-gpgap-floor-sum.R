# Review of lane gpgap: the FLOOR lines of the floor logs, as a table.
#   Rscript dev/review-gpgap-floor-sum.R
fs <- Sys.glob("dev/review-gpgap-out/floor-*.txt")
L <- unlist(lapply(fs, readLines, warn = FALSE))
L <- L[startsWith(L, "FLOOR")]
num <- function(pat) as.numeric(sub(pat, "\\1", L))
g <- function(p) regmatches(L, regexpr(p, L))
tb <- data.frame(
  tag = sub("^FLOOR (\\S+) .*", "\\1", L),
  seed = num(".* seed (\\d+) .*"),
  n = num(".* n (\\d+) A .*"),
  A = num(".* A (\\S+) xs .*"),
  xs = num(".* xs (\\S+) rep .*"),
  rep = num(".* rep (\\d+) ratio .*"),
  ratio = num(".* ratio (\\S+) \\| brms.*"),
  s_b = num(".*brms: sigma (\\S+) sdgp.*"),
  sd_b = num(".*brms: sigma \\S+ sdgp (\\S+) rho.*"),
  rho_b = num(".*brms: sigma \\S+ sdgp \\S+ rho (\\S+) ll.*"),
  s_f = num(".*frmnug: sigma (\\S+) sdgp.*"),
  sd_f = num(".*frmnug: sigma \\S+ sdgp (\\S+) rho.*"),
  rel = num(".*sigma/brms (\\S+), .*"),
  dse = suppressWarnings(num(".*\\|dlog\\|/SE (\\S+), .*")),
  loss = num(".*brms-at-frm loss (\\S+), .*"),
  ll_b = num(".*brms: sigma \\S+ sdgp \\S+ rho \\S+ ll (\\S+) \\| frmnug.*"),
  ll_f = num(".*frmnug: sigma \\S+ sdgp \\S+ rho \\S+ ll (\\S+) \\| sigma/brms.*"),
  sb_over = num(".*sigma_b/sdgp_b (\\S+) \\|.*"),
  frm_over = suppressWarnings(num(".*sigma/sdgp (\\S+) warn.*")),
  frm_s = suppressWarnings(num(".*frm\\(\\): sigma (\\S+) sdgp.*")),
  warned = grepl("nugwarn TRUE", L),
  stringsAsFactors = FALSE)
tb$gap <- tb$ll_b - tb$ll_f
tb <- tb[order(tb$tag, tb$seed, tb$n, -tb$ratio), ]
options(width = 200)
print(transform(tb, s_b = signif(s_b, 4), s_f = signif(s_f, 4), rel = round(rel, 4),
                dse = round(dse, 2), loss = signif(loss, 3), gap = round(gap, 3),
                sb_over = signif(sb_over, 3), frm_over = signif(frm_over, 3),
                frm_s = signif(frm_s, 4))[, c("tag", "seed", "n", "A", "xs", "rep", "ratio",
                 "s_b", "s_f", "frm_s", "sd_b", "sd_f", "rel", "dse", "loss", "gap",
                 "sb_over", "frm_over", "warned")], row.names = FALSE)
cat("\nWarning against effect, all rows with rep = 1:\n")
t1 <- tb[tb$rep == 1, ]
t1$bin <- cut(t1$sb_over, c(0, 3e-3, 5e-3, 7e-3, 1e-2, 1.2e-2, 2e-2, 1), right = FALSE)
print(do.call(rbind, lapply(split(t1, t1$bin), function(z) {
  if (!nrow(z)) return(NULL)
  data.frame(rows = nrow(z), warned = sum(z$warned),
             min_rel = round(min(z$rel), 4), max_rel = round(max(z$rel), 4),
             med_dse = round(median(z$dse, na.rm = TRUE), 2),
             max_dse = round(max(z$dse, na.rm = TRUE), 2),
             max_loss = signif(max(z$loss), 3), max_gap = round(max(z$gap), 3))
})))
