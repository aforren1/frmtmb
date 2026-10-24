# Fix round 1, part B: bare (1 | gr(g)) is (1 | g) as in brms (B1), and
# the bhaz() refusal with its pointer to cox() in one place and its
# gr = sentence (B2, B3), on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
rep_file <- function(p, old, new) {
  p <- file.path(root, p)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
  if (regexpr(old, s, fixed = TRUE) < 0 || n != 1) {
    stop(p, ": count ", n, ": ", substr(old, 1, 70))
  }
  writeLines(sub(old, new, s, fixed = TRUE), p)
}
rep_file("R/parse.R", '      if (length(gvar) != 1 || (has_cov + has_prec) > 1 ||
          (is.null(st) && !has_by && (has_cov + has_prec) != 1) ||
          !all(nms %in% c("", "cov", "prec", "by"))) {
        frm_stop("gr() supports (x | gr(g, cov = A)), ",
                 "(1 | gr(g, prec = Q)), (x | gr(g, by = f)) and ",
                 "(x | gr(g, dist = \\"student\\"))", call. = FALSE)
      }', '      if (length(gvar) != 1 || (has_cov + has_prec) > 1 ||
          !all(nms %in% c("", "cov", "prec", "by"))) {
        frm_stop("gr() supports (x | gr(g)), (x | gr(g, cov = A)), ",
                 "(1 | gr(g, prec = Q)), (x | gr(g, by = f)) and ",
                 "(x | gr(g, dist = \\"student\\"))", call. = FALSE)
      }
      # brms\'s bare gr(g) is the grouping factor g itself (its standata
      # is that of (1 | g)); a multi-membership group keeps its own
      # spelling, which the term above would no longer reach
      if (is.null(st) && !has_by && !has_cov && !has_prec &&
            calls_function(gvar[[1L]], "mm")) {
        frm_stop("gr(mm(...)) is not a spelling of a multi-membership term: ",
                 "write (x | mm(g1, g2)) without gr()", call. = FALSE)
      }')
p <- "R/brms-vocabulary.R"
rep_file(p, '# brms addition terms frmtmb refuses, with the reason and the frmtmb
# spelling where there is one.
brms_aterm_refusals <- c(
  bhaz = paste0(
    "bhaz() sets brms\'s M-spline basis for the baseline hazard of a ",
    "cox() model. frmtmb sets the same basis in the family instead: ",
    "family = cox(df = 5, degree = 3, intercept = TRUE) is brms\'s ",
    "bhaz() default")
)', '# Where the bhaz() refusal sends a user: the family that takes brms\'s
# baseline-hazard basis. The one place to change when cox() moves to an
# extension package (the refusal reads it).
bhaz_family_pointer <- paste0(
  "family = cox(df = 5, degree = 3, intercept = TRUE), which is brms\'s ",
  "bhaz() default")

# brms addition terms frmtmb refuses, with the reason and the frmtmb
# spelling where there is one.
brms_aterm_refusals <- c(
  bhaz = paste0(
    "it sets brms\'s M-spline basis for the baseline hazard of a cox() ",
    "model, and frmtmb sets the same basis in the family: ",
    bhaz_family_pointer, ". A group-specific baseline, bhaz(gr = ), has ",
    "no frmtmb equivalent")
)')
cat("fix-b applied to", root, "\n")
