# Reviewer, lane wt-conditions: printed-output fidelity of the helpers
# against the base functions, one fresh Rscript per case and arm.
# Base arm: base library, STOP/WARN/MSG are stop/warning/message.
# Lane arm: lane library, STOP/WARN/MSG are frm_stop/frm_warning/
# frm_message. Both arms load frmtmb quietly first. Stdout and stderr
# are captured separately and together (2>&1).
#   Rscript dev/conditions-rev-shapes.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
out <- file.path(root, "dev/conditions-rev-log/shapes")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rscript <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
libs <- c(base = "C:/Users/adf44/source/r/rellib-r3",
          lane = "C:/Users/adf44/source/r/conditions-lib")
subs <- list(base = c(STOP = "stop", WARN = "warning", MSG = "message"),
             lane = c(STOP = "frm_stop", WARN = "frm_warning",
                      MSG = "frm_message"))

cases <- list(
  e_direct = 'f <- function(x) STOP("bad x: ", x)\nf(3)',
  e_nocall = 'f <- function() STOP("quiet", call. = FALSE)\nf()',
  e_toplevel = 'STOP("top level")',
  e_try_print = 'f <- function(x) STOP("bad x: ", x)\ntry(f(3))\ncat(geterrmessage())',
  e_try_long = 'a_function_with_a_long_name <- function(argument_one, argument_two) STOP("a message that is long enough to wrap the try prefix")\ntry(a_function_with_a_long_name(argument_one = 1, argument_two = 2))',
  e_promise = 'g <- function(a) a\nf <- function() g(STOP("p"))\nf()',
  e_lapply = 'f <- function() lapply(1, function(i) STOP("l ", i))\nf()',
  e_map = 'f <- function() Map(function(a) STOP("m ", a), 1)\nf()',
  e_vapply = 'f <- function() vapply(1, function(i) STOP("v"), 1)\nf()',
  e_handler = 'f <- function() tryCatch(log("a"), error = function(e) STOP("h: ", conditionMessage(e)))\nf()',
  e_wch = 'f <- function() withCallingHandlers(warning("w0"), warning = function(w) STOP("wch"))\nf()',
  e_onexit = 'f <- function() { on.exit(STOP("oe")); 1 }\nf()',
  e_eval = 'f <- function() eval(quote(STOP("ev")))\nf()',
  e_evalq_env = 'f <- function() { e <- new.env(); evalq(STOP("evq"), e) }\nf()',
  e_docall = 'f <- function(x) STOP("dc")\ndo.call(f, list(1:3))',
  e_docall_big = 'f <- function(x) STOP("dc big")\ntry(do.call(f, list(mtcars)))',
  e_s3 = 'gen <- function(x) UseMethod("gen")\ngen.default <- function(x) STOP("s3")\ngen(1)',
  e_nextmethod = 'gen <- function(x) UseMethod("gen")\ngen.a <- function(x) NextMethod()\ngen.default <- function(x) STOP("s3n")\ngen(structure(1, class = "a"))',
  e_s4 = 'setGeneric("gen4", function(x) standardGeneric("gen4"))\nsetMethod("gen4", "numeric", function(x) STOP("s4"))\ngen4(1)',
  e_s4_arg = 'setGeneric("gen4", function(x) standardGeneric("gen4"))\nsetMethod("gen4", "numeric", function(x) x)\nh <- function() STOP("in arg")\ngen4(h())',
  e_recall = 'f <- function(n) if (n > 0) Recall(n - 1) else STOP("rec")\nf(2)',
  e_multiline = 'f <- function() STOP(sprintf("line one %d\\nline two %s", 1L, "x"), "\\n  tail")\nf()',
  e_gettextf = 'f <- function(n) STOP(gettextf("got %d items", n), domain = NA)\nf(2)',
  e_empty = 'f <- function() STOP()\nf()',
  e_vector = 'f <- function() STOP(c("a", "b"), 1:2)\nf()',
  e_utf8 = 'f <- function() STOP("sigma \\u03c3 \\u2265 0")\nf()',
  e_nested_calls = 'h <- function() STOP("deep")\ng <- function() h()\nf <- function() g()\nf()',
  e_showcalls_off = 'options(showErrorCalls = FALSE)\nh <- function() STOP("deep")\ng <- function() h()\ng()',
  e_warnlen = 'options(warning.length = 100)\nf <- function() STOP(strrep("x", 300))\nf()',
  e_error_opt = 'options(error = function() cat("HOOK", geterrmessage()))\nf <- function() STOP("hooked")\nf()',
  e_catch_error = 'f <- function() STOP("c")\nr <- tryCatch(f(), error = function(e) "caught")\nprint(r)\nprint(inherits(tryCatch(f(), error = identity), "simpleError"))',
  w_deferred = 'f <- function(x) WARN("careful ", x)\nf(1)\ncat("after\\n")',
  w_nocall = 'f <- function(x) WARN("careful ", x, call. = FALSE)\nf(1)\ncat("after\\n")',
  w_toplevel = 'WARN("top warn")\ncat("after\\n")',
  w_warn1 = 'options(warn = 1)\nf <- function(x) WARN("careful ", x)\nf(1)\ncat("after\\n")',
  w_warn2 = 'options(warn = 2)\nf <- function(x) WARN("careful ", x)\nf(1)\ncat("after\\n")',
  w_warn_neg = 'options(warn = -1)\nf <- function(x) WARN("careful ", x)\nf(1)\ncat("after\\n")',
  w_many = 'f <- function(i) WARN("w", i)\ng <- function() for (i in 1:12) f(i)\ng()\nprint(warnings())',
  w_two = 'f <- function(i) WARN("w", i)\ng <- function() { f(1); f(2) }\ng()',
  w_long = 'a_function_with_a_long_name <- function(argument_one, argument_two) WARN("a message that is long enough to wrap the prefix line")\na_function_with_a_long_name(argument_one = 1, argument_two = 2)',
  w_suppress = 'f <- function(x) WARN("careful ", x)\nsuppressWarnings(f(1))\nsuppressWarnings(f(1), classes = "warning")\ncat("after\\n")',
  w_wch_muffle = 'f <- function(x) WARN("careful ", x)\nwithCallingHandlers(f(1), warning = function(w) { cat("got:", conditionMessage(w), deparse(conditionCall(w)), "\\n"); invokeRestart("muffleWarning") })',
  w_trycatch = 'f <- function(x) WARN("careful ", x)\nprint(tryCatch(f(1), warning = function(w) conditionMessage(w)))',
  w_value = 'f <- function(x) WARN("careful ", x)\nv <- suppressWarnings(f(1))\nprint(v)',
  w_immediate = 'f <- function() { WARN("imm", immediate. = TRUE); cat("between\\n") }\nf()',
  w_promise = 'g <- function(a) a\nf <- function() g(WARN("pw"))\nf()',
  w_lapply = 'f <- function() invisible(lapply(1, function(i) WARN("lw")))\nf()',
  w_in_handler_last = 'f <- function(x) WARN("careful ", x)\nf(1)\nprint(names(warnings()))',
  m_basic = 'f <- function() MSG("note ", 1)\nf()\ncat("after\\n")',
  m_nolf = 'f <- function() MSG("no newline", appendLF = FALSE)\nf()\ncat("|after\\n")',
  m_suppress = 'f <- function() MSG("note")\nsuppressMessages(f())\nsuppressMessages(f(), classes = "message")\ncat("after\\n")',
  m_sink = 'f <- function() MSG("note")\nzz <- textConnection("cap", "w")\nsink(zz, type = "message")\nf()\nsink(type = "message")\nclose(zz)\nprint(cap)',
  m_call = 'f <- function() MSG("note")\nprint(conditionCall(tryCatch(f(), message = identity)))',
  m_wch = 'f <- function() MSG("note")\nwithCallingHandlers(f(), message = function(m) { cat("got:", conditionMessage(m)); invokeRestart("muffleMessage") })'
)

run_arm <- function(nm, arm) {
  code <- cases[[nm]]
  for (k in names(subs[[arm]])) code <- gsub(k, subs[[arm]][[k]], code,
                                             fixed = TRUE)
  pre <- sprintf(paste0(".libPaths(c('%s', 'C:/Users/adf44/source/r/",
                        "pinlib', 'C:/Users/adf44/AppData/Local/R/",
                        "win-library/4.6'))\n",
                        "suppressPackageStartupMessages(library(frmtmb))\n"),
                 libs[[arm]])
  f <- file.path(out, sprintf("%s.%s.R", nm, arm))
  writeLines(c(pre, code), f)
  both <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE,
                                   stderr = TRUE))
  so <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE,
                                 stderr = FALSE))
  list(both = both, stdout = so)
}

res <- list()
for (nm in names(cases)) {
  b <- run_arm(nm, "base"); l <- run_arm(nm, "lane")
  same_both <- identical(b$both, l$both)
  same_out <- identical(b$stdout, l$stdout)
  # the only allowed difference: a call that shows the raising call itself
  norm <- function(x) gsub("frm_(stop|warning|message)", "\\1", x)
  same_norm <- identical(norm(b$both), norm(l$both))
  cat(sprintf("%-18s both %-5s stdout %-5s after-name-normalizing %s\n",
              nm, same_both, same_out, same_norm))
  if (!same_both || !same_out) {
    cat("  --- base (2>&1)\n"); cat(paste0("  | ", b$both), sep = "\n")
    cat("  --- lane (2>&1)\n"); cat(paste0("  | ", l$both), sep = "\n")
    if (!same_out) {
      cat("  --- base stdout\n"); cat(paste0("  | ", b$stdout), sep = "\n")
      cat("  --- lane stdout\n"); cat(paste0("  | ", l$stdout), sep = "\n")
    }
  }
  res[[nm]] <- c(both = same_both, stdout = same_out, norm = same_norm)
}
m <- do.call(rbind, res)
cat("\ncases", nrow(m), "identical 2>&1", sum(m[, "both"]),
    "identical stdout", sum(m[, "stdout"]),
    "identical after name normalizing", sum(m[, "norm"]), "\n")
