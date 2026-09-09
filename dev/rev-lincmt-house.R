# Attack 7: housekeeping, and the two guards. A guard is worth nothing
# until the case it guards against has been CONSTRUCTED.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

cat("\n=== A. the registry, exactly one row added ===\n")
ft <- frmtmb::frm_compat_features()
rl <- frmtmb::frm_compat_rules()
print(ft[grepl("frm_ode|frm_lincmt", ft$name), ])
cat("rules naming frm_lincmt():",
    nrow(rl[rl$feature_a == "frm_lincmt()" |
              rl$feature_b == "frm_lincmt()", ]), "\n")
cat("rules naming frm_ode():",
    nrow(rl[rl$feature_a == "frm_ode()" | rl$feature_b == "frm_ode()",
            ]), "\n")

cat("\n=== B. the same registry at the BASE commit, for the diff ===\n")
base <- new.env()
local({
  op <- .libPaths()
  on.exit(.libPaths(op))
}, envir = base)
cat("(run separately below)\n")

cat("\n=== C. the compat assertion, SEEN TO FAIL ===\n")
# the assertion in test-compat.R is nrow(rules naming frm_lincmt) == 0.
# Plant the row it exists to catch and re-evaluate the same predicate.
row <- rl[1L, , drop = FALSE]
row[["feature_a"]] <- "frm_lincmt()"
row[["feature_b"]] <- "frm_sample"
planted <- rbind(rl, row)
pred <- function(tab) nrow(tab[tab$feature_a == "frm_lincmt()" |
                                 tab$feature_b == "frm_lincmt()", ])
cat("  shipped table:", pred(rl), " (the test asserts 0)\n")
cat("  with a planted rule:", pred(planted),
    " -> the assertion would", if (pred(planted) == 0L)
      "STILL PASS: the guard fails open" else "FAIL", "\n")

cat("\n=== D. the hazard scanner, SEEN TO FIRE on lincmt.R ===\n")
ns <- asNamespace("frmtmb.ode")
cat("  shipped:", length(frmtmb::frm_hazard_reads("frmtmb.ode")),
    "hazard reads\n")
orig <- get("lincmt_resp", envir = ns)
body_txt <- paste(deparse(body(orig)), collapse = "\n")
planted_fn <- orig
body(planted_fn) <- as.call(list(
  quote(`{`),
  quote(if (FALSE) frame$data),
  body(orig)))
environment(planted_fn) <- ns
unlockBinding("lincmt_resp", ns)
assign("lincmt_resp", planted_fn, envir = ns)
hits <- frmtmb::frm_hazard_reads("frmtmb.ode")
cat("  with a `$` read planted in lincmt_resp():", length(hits),
    "\n")
if (length(hits)) cat("   ", paste(hits, collapse = "\n    "), "\n")
assign("lincmt_resp", orig, envir = ns)
stopifnot(identical(get("lincmt_resp", envir = ns), orig))
cat("  restored:", length(frmtmb::frm_hazard_reads("frmtmb.ode")),
    "hazard reads\n")

cat("\n=== E. lincmt.R house style ===\n")
p <- file.path("C:/Users/adf44/source/r/frmtmb-wt-lincmt",
               "extensions/frmtmb.ode/R/lincmt.R")
ln <- readLines(p, warn = FALSE)
long <- which(nchar(ln) > 80)
cat("  lines over 80 columns:", length(long),
    if (length(long)) paste(long, collapse = " ") else "", "\n")
cat("  em dashes:", sum(grepl("\u2014", ln)), "\n")
cat("  spaced hyphens standing in for one:",
    sum(grepl("[[:alnum:]] - [[:alnum:]]", ln) &
          !grepl("^#'? *[0-9]|<-|\\+|=", ln)), "\n")
nrd <- grep("@noRd", ln)
bad <- nrd[vapply(nrd, function(i)
  i < length(ln) && grepl("^#'", ln[[i + 1L]]) &&
    nzchar(trimws(sub("^#'", "", ln[[i + 1L]]))), TRUE)]
cat("  @noRd followed by roxygen text:", length(bad), "\n")

cat("\n=== F. test-lincmt.R: any absolute numeric tolerance ===\n")
tp <- file.path("C:/Users/adf44/source/r/frmtmb-wt-lincmt",
                "extensions/frmtmb.ode/tests/testthat/test-lincmt.R")
tl <- readLines(tp, warn = FALSE)
cand <- grep("expect_(lt|gt|equal|lte|gte)", tl, value = TRUE)
cat("  assertions with a bare 1e- constant on the bar side:\n")
for (s in cand) if (grepl("1e-[0-9]+\\s*\\)", s))
  cat("   ", trimws(s), "\n")
cat("  (a 1e- constant used as a SEPARATION rather than a bar is",
    "fine)\n")
cat("  lines over 80 columns in the test file:",
    sum(nchar(tl) > 80), "\n")
