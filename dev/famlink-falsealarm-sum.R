# Join dev/famlink-falsealarm-{base,lane}.tsv and print the generated
# block that dev/famlink-findings.md carries verbatim.
b <- utils::read.delim("dev/famlink-falsealarm-base.tsv", colClasses = "character",
                       quote = "")
l <- utils::read.delim("dev/famlink-falsealarm-lane.tsv", colClasses = "character",
                       quote = "")
m <- merge(b, l, by = c("family", "link"), suffixes = c("_base", "_lane"))
stopifnot(nrow(m) == nrow(b), nrow(m) == nrow(l))
false_alarm <- m$status_base == "fit" & m$status_lane != "fit"
both_fit <- m$status_base == "fit" & m$status_lane == "fit"
same_ll <- both_fit & m$loglik_base == m$loglik_lane
nonconv <- both_fit & m$convergence_lane != "0"
cat("---- GENERATED: dev/famlink-falsealarm-sum.R ----\n")
cat(sprintf("(family, link) pairs brms 2.23.0 accepts, excluding multinomial: %d\n",
            nrow(m)))
cat(sprintf("families: %d\n", length(unique(m$family))))
cat(sprintf("fitted on base: %d; fitted on lane: %d\n",
            sum(m$status_base == "fit"), sum(m$status_lane == "fit")))
cat(sprintf("false alarms (base fitted, lane refused): %d\n", sum(false_alarm)))
cat(sprintf("fitted in both with logLik identical to 10 decimals: %d of %d\n",
            sum(same_ll), sum(both_fit)))
cat(sprintf("fitted in both, lane optimizer code not 0: %d\n", sum(nonconv)))
err <- m[m$status_lane != "fit", c("family", "link", "status_base", "message_lane")]
err$message_lane <- substr(err$message_lane, 1, 60)
cat("pairs the lane does not fit:\n")
print(err, row.names = FALSE)
if (sum(nonconv)) {
  cat("pairs fitted with a nonzero optimizer code (same in base?):\n")
  print(m[nonconv, c("family", "link", "convergence_base", "convergence_lane")],
        row.names = FALSE)
}
cat("---- END GENERATED ----\n")
