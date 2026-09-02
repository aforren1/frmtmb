# Shared library for the brms-vignette port audit.
#
# Why AST rather than regex: the transform must remove *named arguments*
# from specific calls without touching identically named variables or
# arguments of other functions (`control` is a real frm() argument, and
# `seed` appears in set.seed()). Deparsing a modified parse tree also
# guarantees the emitted code is syntactically valid.

BRMS_DOC <- system.file("doc", package = "brms")

# Arguments removed from brm()/brm_multiple()/update() calls. Exactly the
# audit's specified list; anything outside it is left alone so that an
# unhandled MCMC argument shows up as a finding rather than silently
# disappearing.
DROP_ARGS <- c(
  "prior", "priors", "chains", "iter", "warmup", "cores", "backend",
  "threads", "refresh", "seed", "control", "init", "file", "silent",
  "save_pars", "sample_prior", "algorithm", "future", "normalize",
  "stanvars", "stan_funs"
)

RENAME_CALLS <- c(brm = "frm", brm_multiple = "frm_multiple")

# Expressions matching these are diverted to the "brms-only" bucket
# instead of being run: frmtmb has set_prior() but the semantics are
# penalized-likelihood, not a Bayesian prior, so auto-converting would
# fabricate agreement.
PRIOR_RE <- "\\b(set_prior|get_prior|prior_string|validate_prior)\\s*\\(|(^|[^_[:alnum:].])prior\\s*\\("

## ---------------------------------------------------------------- extract

extract_rmd <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^[`~]{3,}\\{[rR][ ,}]", lines)
  out <- list()
  for (s in starts) {
    fence <- sub("^([`~]{3,}).*$", "\\1", lines[s])
    e <- s + which(grepl(paste0("^", fence, "\\s*$"), lines[(s + 1):length(lines)]))[1]
    hdr <- lines[s]
    label <- sub("^[`~]{3,}\\{[rR],?\\s*", "", hdr)
    label <- sub("\\}\\s*$", "", label)
    code <- if (e > s + 1) lines[(s + 1):(e - 1)] else character()
    out[[length(out) + 1]] <- list(header = label, code = code)
  }
  out
}

# The JSS vignettes (overview, multilevel) ship as pre-rendered LaTeX; the
# user-visible code lives in Sinput verbatim blocks with the R> prompt.
extract_ltx <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("\\\\begin\\{Sinput\\}", lines)
  ends <- grep("\\\\end\\{Sinput\\}", lines)
  # Both JSS vignettes wrap printed output in Sinput as well as code, and
  # brms_overview marks only some blocks with the R> prompt. So: a
  # prompted block is code by its prompts; an unprompted block is code
  # only if it parses as R and is not a bare-formula syntax illustration.
  out <- list()
  for (i in seq_along(starts)) {
    s <- starts[i]
    e <- ends[ends > s][1]
    if (is.na(e) || e <= s + 1) next
    raw <- lines[(s + 1):(e - 1)]
    if (any(grepl("^R> ", raw))) {
      raw <- raw[grepl("^R> ", raw) | grepl("^\\+ ", raw)]
      code <- sub("^\\+ ", "  ", sub("^R> ", "", raw))
    } else {
      code <- raw
      p <- tryCatch(parse(text = paste(code, collapse = "\n")),
                    error = function(e) NULL)
      if (is.null(p) || !length(p)) next
      if (all(vapply(p, function(x) is.call(x) && identical(call_name(x), "~"),
                     logical(1)))) next
    }
    out[[length(out) + 1]] <- list(header = "", code = code)
  }
  out
}

extract_vignette <- function(name) {
  rmd <- file.path(BRMS_DOC, paste0(name, ".Rmd"))
  ltx <- file.path(BRMS_DOC, paste0(name, ".ltx"))
  chunks <- if (file.exists(rmd)) extract_rmd(rmd) else extract_ltx(ltx)
  chunks <- Filter(function(z) length(z$code) && any(nzchar(trimws(z$code))), chunks)
  # SETTINGS-* chunks are knitr plumbing, never user-facing code
  chunks <- Filter(function(z) !grepl("SETTINGS", z$header), chunks)
  for (i in seq_along(chunks)) {
    chunks[[i]]$vignette <- name
    chunks[[i]]$idx <- i
  }
  chunks
}

## -------------------------------------------------------------- transform

call_name <- function(e) {
  if (!is.call(e)) return(NA_character_)
  f <- e[[1]]
  if (is.name(f)) return(as.character(f))
  # brms::brm(...)
  if (is.call(f) && identical(as.character(f[[1]]), "::")) return(as.character(f[[3]]))
  NA_character_
}

.dropped <- NULL

walk_expr <- function(e) {
  if (!is.call(e)) return(e)
  nm <- call_name(e)
  if (!is.na(nm) && nm %in% names(RENAME_CALLS)) {
    e[[1]] <- as.name(RENAME_CALLS[[nm]])
    nm <- RENAME_CALLS[[nm]]
  }
  if (!is.na(nm) && nm %in% c("frm", "frm_multiple", "update")) {
    anames <- names(e)
    if (!is.null(anames)) {
      kill <- which(anames %in% DROP_ARGS)
      if (length(kill)) {
        .dropped <<- c(.dropped, anames[kill])
        e <- e[-kill]
      }
    }
  }
  for (i in seq_along(e)) {
    if (!is.null(e[[i]]) && (is.call(e[[i]]))) e[[i]] <- walk_expr(e[[i]])
  }
  e
}

# Returns a list per top-level expression: the transformed source, whether
# it was diverted to the brms-only bucket, and which arguments were removed.
transform_code <- function(code) {
  txt <- paste(code, collapse = "\n")
  exprs <- tryCatch(parse(text = txt, keep.source = TRUE), error = function(e) e)
  if (inherits(exprs, "error")) {
    return(list(list(src = txt, status = "PARSE-ERROR", dropped = character(),
                     msg = conditionMessage(exprs))))
  }
  srcs <- attr(exprs, "srcref")
  lapply(seq_along(exprs), function(i) {
    orig <- if (!is.null(srcs)) paste(as.character(srcs[[i]]), collapse = "\n") else
      paste(deparse(exprs[[i]]), collapse = "\n")
    # Never execute package installation from a vignette
    if (grepl("install_github|install\\.packages|remotes::", orig)) {
      return(list(src = orig, status = "SETUP-SKIP", dropped = character(),
                  msg = "package installation"))
    }
    .dropped <<- character()
    new <- walk_expr(exprs[[i]])
    drops <- .dropped
    changed <- !identical(deparse(new), deparse(exprs[[i]]))
    src <- if (changed) paste(deparse(new), collapse = "\n") else orig
    # Bucketed only if prior code survives the transform: a `prior =`
    # argument inside brm() is removed by the drop list, so those models
    # still run; a standalone set_prior()/get_prior() expression does not.
    if (grepl(PRIOR_RE, src)) {
      return(list(src = src, status = "BRMS-ONLY", dropped = unique(drops),
                  msg = "prior code"))
    }
    list(src = src,
         status = if (changed) "TRANSFORMED" else "VERBATIM",
         dropped = unique(drops), msg = "")
  })
}

## ------------------------------------------------------------ classifying

is_model_call <- function(src) grepl("\\bfrm\\s*\\(|\\bfrm_multiple\\s*\\(", src)

#' Score one merged result row (a row of results-merged.rds).
#'
#' CLEAN is decided by pass 1 alone. An expression that errored under the
#' mechanical transform is not clean, even when the edit that rescued it
#' was made to a different expression: that still costs the porter an
#' edit, so folding it into CLEAN would overstate the headline. Those
#' rows are SPELLING: upstream.
#'
#' @noRd
classify <- function(r) {
  # Never executed by design; not a score
  if (r$status_raw %in% c("BRMS-ONLY", "SETUP-SKIP", "PARSE-ERROR")) {
    return(r$status_raw)
  }
  if (identical(r$status_raw, "OK")) return("CLEAN")
  if (identical(r$status_spell, "OK")) {
    return(if (nzchar(r$patch)) paste0("SPELLING: ", r$patch) else "SPELLING: upstream")
  }
  if (grepl("^object '.*' not found", r$msg_spell)) return("CASCADE")
  "FAIL"
}
