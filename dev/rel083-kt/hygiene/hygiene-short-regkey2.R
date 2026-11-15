# Short check of lane hygiene's fix round 2, F1: can the final
# registration_key() still merge two different checks?
#   Rscript dev/reviews/hygiene-short-regkey2.R
# Each pair: the first check refuses, the second does nothing. "added 1"
# is a merge; "fit: fitted" after a merge means the refusal was lost.
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
src <- new.env()
sys.source("C:/Users/adf44/source/r/frmtmb-wt-rel083/R/registries.R",
           envir = src, keep.source = FALSE)
for (f in c("registration_key", "register_once", "registration_named_env")) {
  cat("INSTALLED", f, "== lane source:",
      identical(deparse(get(f, ns)), deparse(get(f, src))), "\n")
}
ce <- ns$frmtmb_frame_checks
n0 <- length(ce$fns)
set.seed(1)
d <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
pair <- function(label, f1, f2, after = NULL) {
  ce$fns <- ce$fns[seq_len(n0)]
  frmtmb_register_frame_check(f1)
  frmtmb_register_frame_check(f2)
  added <- length(ce$fns) - n0
  if (!is.null(after)) after()
  fit <- tryCatch({ suppressMessages(frm(y ~ x, data = d)); "fitted" },
                  error = function(e) paste("refused:", conditionMessage(e)))
  cat(sprintf("%-58s added %d; fit: %s\n", label, added, fit))
  ce$fns <- ce$fns[seq_len(n0)]
}
in_ns <- function(f) {
  environment(f) <- asNamespace("stats")
  f
}

## a list with an environment deep inside ----------------------------------
mk_list <- in_ns(function(cfg) {
  force(cfg)
  function(spec, frame) {
    if (isTRUE(cfg$a$b$env$refuse)) stop("refused, deep env", call. = FALSE)
  }
})
e1 <- new.env(); e2 <- new.env()
pair("DEEP list(a = list(b = list(env = <env>))), 2 envs",
     mk_list(list(a = list(b = list(env = e1)))),
     mk_list(list(a = list(b = list(env = e2)))),
     after = function() e1$refuse <- TRUE)
pair("DEEP the same list object twice (want 1; 2 is a duplicate)",
     mk_list(list(a = list(b = list(env = e1)))),
     mk_list(list(a = list(b = list(env = e1)))))

## S4, reference class, R6-style --------------------------------------------
setClass("Cfg", representation(refuse = "logical"))
mk_s4 <- in_ns(function(cfg) {
  force(cfg)
  function(spec, frame) if (cfg@refuse) stop("refused, S4", call. = FALSE)
})
pair("S4 slot TRUE vs FALSE",
     mk_s4(new("Cfg", refuse = TRUE)), mk_s4(new("Cfg", refuse = FALSE)))
pair("S4 equal slots twice (want 1)",
     mk_s4(new("Cfg", refuse = FALSE)), mk_s4(new("Cfg", refuse = FALSE)))
Rc <- setRefClass("Rc", fields = list(refuse = "logical"))
mk_rc <- in_ns(function(cfg) {
  force(cfg)
  function(spec, frame) if (isTRUE(cfg$refuse)) stop("refused, RC",
                                                     call. = FALSE)
})
r1 <- Rc$new(refuse = FALSE); r2 <- Rc$new(refuse = FALSE)
pair("RC two objects, equal fields, first set to refuse after",
     mk_rc(r1), mk_rc(r2), after = function() r1$refuse <- TRUE)
r6 <- function() structure(new.env(), class = c("Cfg6", "R6"))
q1 <- r6(); q2 <- r6()
pair("R6-style two objects, equal, first set to refuse after",
     mk_rc(q1), mk_rc(q2), after = function() q1$refuse <- TRUE)

## a chain through unnamed environments above the first frame --------------
chain <- function(mode) {
  top <- new.env(parent = globalenv())
  top$mode <- mode
  mid <- new.env(parent = top)
  own <- new.env(parent = mid)
  f <- function(spec, frame) if (mode == "refuse") stop("refused, chain",
                                                         call. = FALSE)
  environment(f) <- own
  f
}
pair("CHAIN value two unnamed frames up differs", chain("refuse"),
     chain("log"))
c1 <- chain("log"); c2 <- chain("log")
pair("CHAIN equal values, first chain turned to refuse after",
     c1, c2, after = function() assign("mode", "refuse",
                                       envir = parent.env(parent.env(
                                         environment(c1)))))
mk_state <- in_ns(function() {
  state <- FALSE
  function(spec, frame) if (state) stop("refused, state", call. = FALSE)
})
s1 <- mk_state(); s2 <- mk_state()
pair("STATE own frame, first set to refuse after",
     s1, s2, after = function() assign("state", TRUE, environment(s1)))

## an environment with a name that is not a namespace ----------------------
named <- function(mode, nm) {
  e <- new.env(parent = globalenv())
  attr(e, "name") <- nm
  e$mode <- mode
  f <- function(spec, frame) if (mode == "refuse") stop("refused, named",
                                                         call. = FALSE)
  environment(f) <- e
  f
}
pair("NAMED own env named 'cfg', values differ",
     named("refuse", "cfg"), named("log", "cfg"))
pair("NAMED own env named 'cfgA' vs 'cfgB', values differ",
     named("refuse", "cfgA"), named("log", "cfgB"))
named_parent <- function(mode) {
  e <- named(mode, "cfg")
  f <- function(spec, frame) if (mode == "refuse") stop("refused, nparent",
                                                         call. = FALSE)
  environment(f) <- new.env(parent = environment(e))
  f
}
pair("NAMED parent named 'cfg', values differ",
     named_parent("refuse"), named_parent("log"))

## values serialize() might write alike -------------------------------------
mk_val <- in_ns(function(x) {
  force(x)
  function(spec, frame) {
    if (is.double(x) && length(x) == 1L && !is.na(x) && 1 / x < 0) {
      stop("refused, -0", call. = FALSE)
    }
  }
})
bytes <- function(x) serialize(x, NULL)
cat("VAL -0 vs 0 bytes equal:", identical(bytes(-0), bytes(0)), "\n")
pair("VAL -0 vs 0", mk_val(-0), mk_val(0))
cat("VAL NA_real_ vs NaN bytes equal:", identical(bytes(NA_real_), bytes(NaN)),
    "\n")
pair("VAL NA_real_ vs NaN", mk_val(NA_real_), mk_val(NaN))
nan2 <- -(0 / 0)
cat("VAL NaN vs -(0/0) bytes equal:", identical(bytes(NaN), bytes(nan2)),
    "; identical():", identical(NaN, nan2), "\n")
pair("VAL NaN vs -(0/0) (want: not merged or harmless)",
     mk_val(NaN), mk_val(nan2))
pair("VAL 1:3 vs c(1L, 2L, 3L) (equal; 2 is a duplicate)",
     mk_val(1:3), mk_val(c(1L, 2L, 3L)))
u <- "\u00e9"; l <- iconv(u, "UTF-8", "latin1")
pair("VAL UTF-8 vs latin1 of one string (equal; 2 is a duplicate)",
     mk_val(u), mk_val(l))
pair("VAL 1L vs 1 (differ in type)", mk_val(1L), mk_val(1))

## the .onLoad() idioms -----------------------------------------------------
onload_named <- in_ns(function(libname, pkgname) {
  check <- function(spec, frame) invisible(NULL)
  frmtmb::frmtmb_register_frame_check(check)
})
ce$fns <- ce$fns[seq_len(n0)]
onload_named("L", "p"); onload_named("L", "p")
cat("ONLOAD `check <- function...; register(check)` twice: added",
    length(ce$fns) - n0, "(1 = merged, 2 = duplicate)\n")
onload_anon <- in_ns(function(libname, pkgname) {
  frmtmb::frmtmb_register_frame_check(function(spec, frame) invisible(NULL))
})
ce$fns <- ce$fns[seq_len(n0)]
onload_anon("L", "p"); onload_anon("L", "p")
cat("ONLOAD anonymous closure twice: added", length(ce$fns) - n0, "\n")
onload_fml <- in_ns(function(libname, pkgname) {
  f <- y ~ x
  frmtmb::frmtmb_register_frame_check(function(spec, frame) invisible(f))
})
ce$fns <- ce$fns[seq_len(n0)]
onload_fml("L", "p"); onload_fml("L", "p")
cat("ONLOAD closure beside a formula twice: added", length(ce$fns) - n0,
    "\n")
ce$fns <- ce$fns[seq_len(n0)]
