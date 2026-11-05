# Generate the measured blocks of dev/nlfix-findings.md from the local logs,
# so that no count or constant is typed. Rscript dev/nlfix-report.R
# Replaces each <!-- nlfix:<name> begin --> ... end block.
num <- function(x) if (length(x)) format(x, digits = 3) else "NA"
read_l1 <- function(dirs) {
  rows <- list()
  for (d in dirs) {
    for (f in list.files(d, "[.]txt$", full.names = TRUE)) {
      ln <- readLines(f, warn = FALSE)
      g <- function(pat) {
        m <- grep(paste0("^ *[$] ", pat, " *:"), ln, value = TRUE)
        if (!length(m)) return(NA_character_)
        trimws(sub("^[^:]*:", "", m[1L]))
      }
      nm <- sub("[.]txt$", "", basename(f))
      rows[[nm]] <- data.frame(
        spec = nm,
        bf = sub("^bf: ", "", grep("^bf: ", ln, value = TRUE)[1L]),
        family = sub("^family: ", "", grep("^family: ", ln, value = TRUE)[1L]),
        const = sub("^num ", "", g("l1_const")),
        grad = sub("^num ", "", g("l1_grad")),
        extra = sub("^num \\[1:2\\] ", "", g("l1_extra_const")),
        l1_ok = sub("^logi ", "", g("l1_ok")),
        verdict = gsub("^chr |\"", "", g("verdict")),
        l1_error = substr(gsub("^chr |\"", "", g("l1_error")), 1, 60),
        brms_error = substr(gsub("^chr |\"", "", g("brms_error")), 1, 70),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
tab_md <- function(df) {
  df[] <- lapply(df, function(v) {
    v <- ifelse(is.na(v), "", v)
    gsub("|", "\\|", v, fixed = TRUE)
  })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
l1 <- read_l1(c("dev/nlfix-l1-out", "dev/nlfix-l3-out"))
l1 <- l1[!grepl("^ctl_|_m$|_fit$|ac_ar2_cov|^mm_dirichlet2$", l1$spec), ]
# the autocor terms' L1 came with their methods run (the same fit, the same
# program)
acl <- read_l1("dev/nlfix-l2-out")
acl <- acl[grepl("^ac_", acl$spec), ]
acl$spec <- sub("_m$", "", acl$spec)
l1 <- rbind(l1, acl)
l1 <- l1[!duplicated(l1$spec), ]
l1 <- l1[order(l1$spec), c("spec", "family", "const", "grad", "extra", "l1_ok",
                           "verdict", "l1_error")]
blocks <- list()
blocks[["l1"]] <- c(
  sprintf("%d specs (`dev/nlfix-l1-specs.R`, through `dev/nlfix-l1.R`, the equiv instrument unchanged).",
          nrow(l1)), "", tab_md(l1))
m <- read_l1("dev/nlfix-l2-out")
ms <- list()
for (f in list.files(c("dev/nlfix-l2-out", "dev/nlfix-l3-out"), "_(m|fit)[.]txt$",
                     full.names = TRUE)) {
  ln <- readLines(f, warn = FALSE)
  at <- grep("^ +[.][.][$] (fitted|simulate|predict) *:", ln)
  meth <- trimws(sub(":.*", "", sub("^ +[.][.][$] ", "", ln[at])))
  val <- vapply(at, function(i) {
    w <- ln[i:min(i + 4L, length(ln))]
    o <- sub(".*logi ", "", grep("[$] ok", w, value = TRUE)[1L])
    me <- sub(".*num ", "", grep("[$] metric", w, value = TRUE)[1L])
    paste0(o, " (", me, ")")
  }, "")
  r <- as.list(setNames(val, meth))
  ms[[f]] <- data.frame(spec = sub("[.]txt$", "", basename(f)),
                        fitted = r$fitted %||% "",
                        simulate = r$simulate %||% "",
                        predict = r$predict %||% "")
}
ms <- do.call(rbind, ms)
ms <- ms[order(ms$spec), ]
rownames(ms) <- NULL
blocks[["l2"]] <- c(sprintf("%d specs. Each cell is ok (metric).", nrow(ms)), "",
                    tab_md(ms))
gen_lines <- function(f, pat = NULL) {
  if (!file.exists(f)) return(paste("missing:", f))
  x <- readLines(f, warn = FALSE)
  if (!is.null(pat)) x <- grep(pat, x, value = TRUE)
  x
}
fence <- function(x) c("```", substr(x, 1, 180), "```")
blocks[["d1check"]] <- fence(gen_lines("dev/nlfix-d1check.log",
                                       "linear -?[0-9]|thres[(]3"))
blocks[["starts"]] <- fence(gen_lines("dev/nlfix-startcheck.log",
                                      "linear +-?[0-9NA]|^sum|^[-0-9]"))
blocks[["ext"]] <- fence(c(gen_lines("dev/nlfix-extfam.log", "linear +-?[0-9]"),
                           gen_lines("dev/nlfix-asym.log",
                                     "linear -?[0-9]|^cox|^rp"),
                           gen_lines("dev/nlfix-sar.log", "logLik|fitted|stancode")))
blocks[["arnl"]] <- fence(gen_lines("dev/nlfix-arnl2.log", ":"))
blocks[["cost2"]] <- fence(gen_lines("dev/nlfix-cost2.log", "guard|saddle fit"))
blocks[["stanc"]] <- fence(gen_lines("dev/nlfix-stanc.log", "^in Stan|^not in Stan"))
blocks[["stanfuns"]] <- fence(gen_lines("dev/nlfix-stanfuns.log"))
sw <- "dev/nlfix-equiv-out/nl-sweep.tsv"
if (file.exists(sw)) {
  s <- read.delim(sw)
  blocks[["sweep"]] <- c(
    tab_md(as.data.frame(table(outcome = s$outcome), responseName = "families")),
    "", tab_md(s[, intersect(c("family", "status_now", "frm_linear", "frm_nl",
                               "frm_intercept", "outcome", "proposal",
                               "default_linear", "default_nl_default_start"),
                             names(s))]))
}
ef <- "dev/nlfix-equiv-compare.txt"
if (file.exists(ef)) blocks[["verdicts"]] <- readLines(ef, warn = FALSE)
for (tag in c("r16", "core", "gated", "ext", "sess", "check")) {
  f <- paste0("dev/nlfix-summary-", tag, ".txt")
  if (file.exists(f)) blocks[[paste0("suite_", tag)]] <- fence(readLines(f, warn = FALSE))
}
sd <- list.files("dev/nlfix-saddle-out", "[.]txt$", full.names = TRUE)
if (length(sd)) {
  t1 <- do.call(rbind, lapply(sd, utils::read.delim))
  ids <- readLines("dev/nlfix-saddle-ids.txt", warn = FALSE)
  miss <- setdiff(ids, t1$id)
  fm <- function(v) format(round(v, 4), nsmall = 4, trim = TRUE)
  t2 <- data.frame(family = t1$family, body = t1$body, default = fm(t1$default),
                   started = fm(t1$started), no_body = fm(t1$no_body),
                   brms_best = fm(t1$brms_best),
                   brms_minus_frm = format(signif(t1$gap_brms, 2), trim = TRUE),
                   restarted = t1$restarted, conv = t1$conv,
                   err = substr(t1$err, 1, 60))
  good <- sum(!is.na(t1$gap_brms) & t1$gap_brms < 1e-3 * abs(t1$brms_best))
  blocks[["saddle"]] <- c(
    sprintf(paste0("%d of %d cases have a row. frm() from the default start is at or above ",
                   "brms's best of 6 on %d; `restarted` is %d. No row: %s."),
            nrow(t1), length(ids), good, sum(t1$restarted),
            if (length(miss)) paste(miss, collapse = ", ") else "none"), "", tab_md(t2))
}
if (file.exists("dev/nlfix-a8cmp.log")) blocks[["a8"]] <- fence(readLines("dev/nlfix-a8cmp.log", warn = FALSE))
if (file.exists("dev/nlfix-rev-a4c.log")) blocks[["a4c"]] <- fence(gen_lines("dev/nlfix-rev-a4c.log", "slope|max row"))
if (file.exists("dev/nlfix-rev-a2.log")) blocks[["a2"]] <- fence(gen_lines("dev/nlfix-rev-a2.log", "ll=|ERROR"))
md <- readLines("dev/nlfix-findings.md", warn = FALSE)
for (nm in names(blocks)) {
  a <- which(md == sprintf("<!-- nlfix:%s begin -->", nm))
  b <- which(md == sprintf("<!-- nlfix:%s end -->", nm))
  if (length(a) != 1L || length(b) != 1L) {
    message("no block ", nm)
    next
  }
  md <- c(md[seq_len(a)], blocks[[nm]], md[b:length(md)])
}
writeLines(md, "dev/nlfix-findings.md")
