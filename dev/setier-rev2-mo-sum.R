# Reviewer of lane setier: the mo() study on the trial merge (optima's
# simplex chart), counts only (the script's held-coordinate reference
# assumes the softmax chart, so its SE comparison is not read here).
rd <- function(arm) do.call(rbind, lapply(1:4, function(c)
  read.delim(sprintf("dev/setier-rev2-log/mo/mo-%s-%d.tsv", arm, c),
             stringsAsFactors = FALSE, quote = "")))
for (arm in c("lane")) {
  d <- rd(arm)
  lost <- d$n_lost > 0 | d$se_finite < d$n_par
  cat(sprintf(paste0("%-5s rows %d | lost-SE fits %d | SE warning %d | ",
                     "code!=0 %d | silent lost %d | warn without loss %d | ",
                     "ce print fail %d | lost reasons: %s\n"), arm, nrow(d),
              sum(lost), sum(d$se_warn & lost), sum(d$code != 0 & lost),
              sum(lost & d$n_warn_fit == 0), sum(d$se_warn & !lost),
              sum(!d$ce_ok),
              paste(names(table(unlist(regmatches(d$lost,
                gregexpr("=[a-z]+", d$lost))))),
                table(unlist(regmatches(d$lost, gregexpr("=[a-z]+",
                                                           d$lost)))),
                collapse = " ")))
}
