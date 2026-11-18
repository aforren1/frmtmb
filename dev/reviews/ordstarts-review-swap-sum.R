# Review of lane ordstarts: summary of dev/reviews/ordstarts-review-swap-*
# against the lane's logs (dev/ordstarts-out/main-*.tsv).
#   Rscript dev/reviews/ordstarts-review-swap-sum.R
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "")
X <- do.call(rbind, lapply(Sys.glob("dev/ordstarts-out/main-*-int-*.tsv"),
                           rd))
X$key <- paste(X$family, X$n, X$seed)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
S <- do.call(rbind, lapply(Sys.glob("dev/reviews/ordstarts-review-swap-c*.tsv"),
                           rd))
S$key <- paste(S$family, S$n, S$seed)
S$llA_lane <- A$ll[match(S$key, A$key)]
S$llB <- B$ll[match(S$key, B$key)]
S$elA_lane <- A$elapsed[match(S$key, A$key)]
S$elB <- B$elapsed[match(S$key, B$key)]
S$tseB <- with(B[match(S$key, B$key), ],
               suppressWarnings(pmax(tse1, tse2, na.rm = TRUE)))
S$degB <- B$degen_warn[match(S$key, B$key)]
cat("rows", nrow(S), "; A reproduces the lane's (|d| < 1e-6):",
    sum(abs(S$llA - S$llA_lane) < 1e-6, na.rm = TRUE), "of",
    sum(!is.na(S$llA)), "\n")
up <- S$llB - S$llA_lane > 1e-3
sbest <- pmax(S$llA, ifelse(S$codeS %in% 0 & S$clsS %in% "proper" &
                              !(S$tseS > 100), S$llS, -Inf), na.rm = TRUE)
S$reach <- sbest >= S$llB - 1e-3
S$above <- sbest > S$llB + 1e-3
S$swap_up <- sbest > S$llA + 1e-3
by <- paste(S$family, S$n)
cat("\nper design: B up moves; of those the swap (proper, code 0,",
    "threshold SE <= 100) reaches B; swap above A anywhere; swap above",
    "B\n")
print(data.frame(
  seeds = tapply(!is.na(S$llA), by, sum),
  B_up = tapply(up, by, sum, na.rm = TRUE),
  swap_reaches_B_on_up = tapply(up & S$reach, by, sum, na.rm = TRUE),
  swap_up = tapply(S$swap_up, by, sum, na.rm = TRUE),
  swap_above_B = tapply(S$above, by, sum, na.rm = TRUE),
  swap_up_where_B_not = tapply(S$swap_up & !up, by, sum, na.rm = TRUE)))
cat("\nB up moves the swap misses:\n")
print(S[up & !S$reach, c("key", "llA", "llB", "llS", "codeS", "clsS",
                         "tseS", "degB")], row.names = FALSE)
cat("\nswap above B by > 1e-3:\n")
print(S[S$above %in% TRUE, c("key", "llA", "llB", "llS", "clsS", "tseS")],
      row.names = FALSE)
cat("\ncost: swap refit elapsed median", median(S$elS, na.rm = TRUE),
    "s; A median", median(S$elA, na.rm = TRUE), "s; B median (lane)",
    median(S$elB, na.rm = TRUE), "s\n")
cat("\nsummary(A): warnings on", sum(S$summwarnA > 0, na.rm = TRUE), "of",
    sum(!is.na(S$summwarnA)), "; messages on", sum(S$summmsgA > 0,
                                                  na.rm = TRUE), "\n")
cat("summary(A) warnings on up-move seeds with A silent at fit:",
    sum(up & S$fitwarnA == 0 & S$summwarnA > 0, na.rm = TRUE), "\n")
w <- S[(S$summwarnA > 0) %in% TRUE, c("key", "summwarnsA")]
if (nrow(w)) print(w, row.names = FALSE)
