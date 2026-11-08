# Lane latperf: which pairs of the equivalence record name an autocorrelation
# term, and what their L2 level reads.
r <- utils::read.delim("dev/equiv-data/results.tsv")
pat <- "^(ar|ma|arma|cosy|unstr)[(]"
k <- grepl(pat, r$feature_a) | grepl(pat, r$feature_b)
print(r[k, c("id", "feature_a", "feature_b", "status", "verdict", "level",
             "l2_what", "l2_metric")], width = 200)
