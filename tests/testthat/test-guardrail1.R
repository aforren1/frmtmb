# Guardrail 1 of dev/core-boundary.md: no silent gaps. Every brms name
# (brms 2.23.0's families, addition terms and formula specials,
# frmtmb:::brms_vocabulary) resolves in core, resolves in an installed
# extension, or is refused by name.
#
# The classification runs in a CHILD process holding core alone, so it
# does not depend on which extensions an earlier file of this session
# loaded, and it is the state a user without the extension is in. Then,
# for each extension a refusal names, a child with that package loaded
# runs the same calls and must not refuse them: a refusal that fires on
# a correct call once the package is there would be worse than none. An
# extension that is not installed skips its arm, so the file runs in
# core's suite alone.
#
# The test is generic: it reads the refusals the code gives, not a list
# of which names live where, so a move that adds a home needs no edit
# here.

# One call or more per brms name, in brms's plainest spelling first and
# then the others brms documents. dry_run = "spec" parses the model and
# stops, so no data or matrix need exist. `ext` is what the call must
# give once the package a refusal names is loaded: "ok", or "refused"
# for a spelling frmtmb refuses by name there (car(M), brms's
# deprecated one location per row).
guard_probes <- function() {
  sp <- function(fo, fam = "gaussian()") {
    sprintf("frm(bf(%s), family = %s, data = NULL, dry_run = 'spec')",
            fo, fam)
  }
  fams <- frmtmb:::brms_vocabulary$families
  out <- rbind(
    data.frame(kind = "family", name = fams,
               code = sprintf("brmsfamily('%s')", fams), ext = "ok"),
    data.frame(kind = "family_string", name = fams,
               code = sprintf(paste0("frm(bf(y ~ x), family = '%s', ",
                                     "data = NULL, dry_run = 'spec')"),
                              fams), ext = "ok"))
  at <- rbind(
    c("bhaz", sp("y | bhaz() ~ x")),
    c("bhaz", sp("y | bhaz(df = 5) ~ x")),
    c("cat", sp("y | cat(k) ~ x", "cumulative()")),
    c("cens", sp("y | cens(ce) ~ x")),
    c("cens", sp("y | cens(ce, y2) ~ x")),
    c("dec", sp("y | dec(d) ~ x")),
    c("index", sp("y | index(i) ~ x")),
    c("mi", sp("y | mi() ~ x")),
    c("mi", sp("y | mi(sdy) ~ x")),
    c("rate", sp("y | rate(r) ~ x", "poisson()")),
    c("se", sp("y | se(s) ~ x")),
    c("se", sp("y | se(s, sigma = TRUE) ~ x")),
    c("subset", sp("y | subset(s) ~ x")),
    c("thres", sp("y | thres(k) ~ x", "cumulative()")),
    c("trials", sp("y | trials(n) ~ x", "binomial()")),
    c("trunc", sp("y | trunc(lb = 0) ~ x")),
    c("trunc", sp("y | trunc(ub = 9) ~ x")),
    c("vint", sp("y | vint(v) ~ x")),
    c("vreal", sp("y | vreal(v) ~ x")),
    c("weights", sp("y | weights(w) ~ x")))
  sx <- rbind(
    c("s", sp("y ~ s(x)")), c("s", sp("y ~ s(x, z)")),
    c("t2", sp("y ~ t2(x, z)")),
    c("gp", sp("y ~ gp(x)")), c("gp", sp("y ~ gp(x, k = 10)")),
    c("gp", sp("y ~ gp(x, z)")),
    c("mo", sp("y ~ mo(o)")), c("me", sp("y ~ me(x, sdx)")),
    c("mi", paste0("frm(bf(y ~ mi(x)) + bf(x | mi() ~ z) + ",
                   "set_rescor(FALSE), data = NULL, dry_run = 'spec')")),
    c("cs", sp("y ~ cs(x)", "sratio()")),
    c("cse", sp("y ~ cse(x)", "sratio()")),
    c("car", sp("y ~ x + car(M)")),
    c("car", sp("y ~ x + car(M, gr = g)")),
    c("car", sp("y ~ x + car(M, gr = g, type = 'icar')")),
    c("sar", sp("y ~ x + sar(M)")), c("fcor", sp("y ~ x + fcor(M)")),
    c("ar", sp("y ~ x + ar()")), c("ar", sp("y ~ x + ar(time = t, gr = g)")),
    c("ma", sp("y ~ x + ma()")), c("arma", sp("y ~ x + arma()")),
    c("cosy", sp("y ~ x + cosy()")),
    c("unstr", sp("y ~ x + unstr(t, g)")),
    c("gr", sp("y ~ x + (1 | gr(g))")),
    c("gr", sp("y ~ x + (1 | gr(g, by = f))")),
    c("gr", sp("y ~ x + (1 + x | gr(g, cor = FALSE))")),
    c("gr", sp("y ~ x + (1 | gr(g, id = 'a'))")),
    c("mm", sp("y ~ x + (1 | mm(g1, g2))")),
    c("mmc", sp("y ~ x + (mmc(x1, x2) | mm(g1, g2))")))
  ext <- ifelse(sx[, 2] == sp("y ~ x + car(M)"), "refused", "ok")
  rbind(out,
        data.frame(kind = "aterm", name = at[, 1], code = at[, 2],
                   ext = "ok"),
        data.frame(kind = "special", name = sx[, 1], code = sx[, 2],
                   ext = ext))
}

# Run every probe in a child with frmtmb and `pkgs` loaded; one line per
# probe, the outcome and the message.
guard_run <- function(probes, pkgs = character()) {
  f <- tempfile(fileext = ".R")
  pf <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f, pf)), add = TRUE)
  saveRDS(probes, pf)
  code <- c(sprintf(".libPaths(%s)", paste0(deparse(.libPaths()),
                                           collapse = "")),
            "suppressMessages(library(frmtmb))",
            sprintf("for (p in %s) suppressMessages(library(p, %s))",
                    paste0(deparse(pkgs), collapse = ""),
                    "character.only = TRUE"),
            sprintf("pr <- readRDS(%s)", deparse(pf)),
            "for (i in seq_len(nrow(pr))) {",
            "  r <- tryCatch({",
            "    suppressWarnings(suppressMessages(",
            "      eval(str2lang(pr$code[i]))))",
            "    'OK'",
            "  }, error = function(e) {",
            "    paste('ERR', gsub('[\\r\\n]+', ' ', conditionMessage(e)))",
            "  })",
            "  cat('PROBE', i, r, '\\n', sep = '\\t')",
            "}",
            "cat('CHILDOK\\n')")
  writeLines(code, f)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", shQuote(f)),
                                  stdout = TRUE, stderr = TRUE))
  hit <- grep("^PROBE\t", out, value = TRUE)
  parts <- strsplit(hit, "\t", fixed = TRUE)
  res <- rep(NA_character_, nrow(probes))
  res[as.integer(vapply(parts, `[[`, "", 2L))] <-
    trimws(vapply(parts, `[[`, "", 3L))
  attr(res, "ok") <- any(out == "CHILDOK")
  attr(res, "out") <- out
  res
}

# What one outcome says about the name: "core" (it resolved), the
# extension a refusal names, "refused" (a refusal that names the
# spelling and says why), or "generic" (a catch-all that names no home
# and no reason: the gap this guardrail exists to close).
guard_classify <- function(res, name) {
  if (identical(res, "OK")) return("core")
  m <- regmatches(res, regexec(perl = TRUE,
    "is (?:in|from) the extension package ([A-Za-z0-9.]+), which is not",
    res))[[1L]]
  if (length(m) == 2L) return(paste0("ext:", m[[2L]]))
  generic <- c("is not a supported family\\. Supported families are",
               "is not supported \\(supported:",
               "could not find function", "object '.*' not found",
               "is not a covariance structure")
  if (any(vapply(generic, grepl, NA, x = res))) return("generic")
  if (!grepl(name, res, fixed = TRUE)) return("generic")
  "refused"
}

test_that("frmtmb's brms vocabulary is brms 2.23.0's", {
  skip_if_not_installed("brms", "2.23.0")
  bns <- asNamespace("brms")
  fams <- sub("^[.]family_", "", ls(bns, all.names = TRUE,
                                     pattern = "^[.]family_"))
  # .family_info() is brms's lookup helper, not a family
  expect_setequal(frmtmb:::brms_vocabulary$families,
                  setdiff(fams, "info"))
  ats <- sub("^resp_", "", grep("^resp_", getNamespaceExports("brms"),
                                value = TRUE))
  expect_setequal(frmtmb:::brms_vocabulary$aterms, ats)
  expect_true(all(frmtmb:::brms_vocabulary$specials %in%
                    getNamespaceExports("brms")))
})

test_that("every brms name resolves in core or an extension, or is refused
           by name", {
  skip_on_cran()
  pr <- guard_probes()
  voc <- frmtmb:::brms_vocabulary
  # the probes cover the whole vocabulary, so a name added to it is
  # checked rather than skipped
  expect_setequal(pr$name[pr$kind == "family"], voc$families)
  expect_setequal(pr$name[pr$kind == "aterm"], voc$aterms)
  expect_setequal(pr$name[pr$kind == "special"], voc$specials)
  res <- guard_run(pr)
  expect_true(attr(res, "ok"), info = paste(attr(res, "out"),
                                            collapse = "\n"))
  cls <- mapply(guard_classify, res, pr$name)
  gaps <- which(cls == "generic" | is.na(res))
  expect_identical(length(gaps), 0L,
                   info = paste(pr$kind[gaps], pr$name[gaps], res[gaps],
                                sep = ": ", collapse = "\n"))
  # the three outcomes all occur, so none of the arms below is vacuous
  expect_true(any(cls == "core"))
  expect_true(any(cls == "refused"))
  expect_true(any(startsWith(cls, "ext:")))
  # every refusal is one the code declares with its reason, so a probe
  # that a slip of its own spelling refuses cannot pass for a refusal.
  # A spelling frmtmb refuses by design wherever the term lives (car(M),
  # brms's deprecated default) is refused and is not counted there
  expect_true(all(cls[pr$ext == "refused"] == "refused" |
                    startsWith(cls[pr$ext == "refused"], "ext:")))
  refused <- unique(pr$name[cls == "refused" & pr$ext == "ok"])
  declared <- c(names(frmtmb:::brms_family_refusals),
                names(frmtmb:::brms_aterm_refusals),
                names(frmtmb:::brms_special_refusals))
  expect_setequal(refused, declared)
  # the classifier tells a generic refusal from a named one
  expect_identical(guard_classify(
    "notafam is not a supported family. Supported families are: gaussian",
    "notafam"), "generic")
  expect_identical(guard_classify(
    "Addition term `zz()` is not supported (supported: weights())", "zz"),
    "generic")

  # the extension arm: each package a refusal names resolves the same
  # calls once it is loaded, except a spelling it refuses by name there
  homes <- unique(sub("^ext:", "", cls[startsWith(cls, "ext:")]))
  absent <- character()
  for (pkg in homes) {
    if (!nzchar(system.file(package = pkg))) {
      absent <- c(absent, pkg)
      next
    }
    at <- which(cls == paste0("ext:", pkg))
    res2 <- guard_run(pr[at, , drop = FALSE], pkg)
    expect_true(attr(res2, "ok"), info = pkg)
    cls2 <- mapply(guard_classify, res2, pr$name[at])
    want <- ifelse(pr$ext[at] == "ok", "core", "refused")
    bad <- which(cls2 != want | is.na(res2))
    expect_identical(length(bad), 0L,
                     info = paste(pkg, pr$name[at][bad], res2[bad],
                                  sep = ": ", collapse = "\n"))
  }
  # a core-only run says which arms it could not run
  if (length(absent)) {
    skip(paste("the extension arm needs these packages installed:",
               paste(absent, collapse = ", ")))
  }
})

test_that("an unloaded extension's addition term is refused by name", {
  # dec() is brms's and frmtmb.eam's. This session may have loaded
  # frmtmb.eam, so the refusal is read off the table that raises it
  expect_identical(unname(frmtmb:::extension_aterm_homes[["dec"]]),
                   "frmtmb.eam")
  expect_error(frmtmb:::check_unsupported_aterm("dec"),
               "dec() addition term is in the extension package frmtmb.eam",
               fixed = TRUE)
  expect_error(frmtmb:::check_unsupported_aterm("bhaz"),
               "family = cox(df = 5", fixed = TRUE)
  # and a term no table names falls through to the general refusal
  expect_null(frmtmb:::check_unsupported_aterm("zz"))
})

test_that("a column named like an extension's special is not the special", {
  # the refusal reads CALL positions; a variable that shares the name is
  # data
  local_frmtmb_registries()
  homes <- frmtmb:::extension_special_homes
  expect_null(frmtmb:::check_extension_specials(quote(car + x)))
  if (length(homes)) {
    nm <- names(homes)[[1L]]
    if (!nm %in% cs_names()) {
      expect_error(frmtmb:::check_extension_specials(
        str2lang(paste0("x + ", nm, "(M, gr = g)"))),
        paste0("is in the extension package ", homes[[1L]]), fixed = TRUE)
    }
  }
})

test_that("bare gr(g) is the grouping factor g, as in brms", {
  # brms's standata for (1 | gr(g)) is that of (1 | g); frmtmb 0.73.0
  # refused it
  set.seed(3)
  d <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
  d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + rnorm(60, 0, 0.5)
  f1 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  f2 <- frm(bf(y ~ x + (1 | gr(g))), family = gaussian(), data = d)
  f3 <- frm(bf(y ~ x + (1 + x | gr(g))), family = gaussian(), data = d)
  f4 <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)
  expect_identical(as.numeric(logLik(f2)), as.numeric(logLik(f1)))
  expect_identical(as.numeric(logLik(f3)), as.numeric(logLik(f4)))
  expect_identical(names(VarCorr(f2)), names(VarCorr(f1)))
  expect_identical(variables(f2), variables(f1))
  # the absent case: the gr() spellings with arguments keep their meaning
  expect_error(frm(bf(y ~ x + (1 | gr(mm(g, g)))), family = gaussian(),
                   data = d), "without gr()", fixed = TRUE)
  expect_error(frm(bf(y ~ x + (1 | gr(g, nope = 1))), family = gaussian(),
                   data = d), "gr() supports", fixed = TRUE)
})

test_that("the bhaz() refusal names what has no equivalent", {
  expect_error(frmtmb:::check_unsupported_aterm("bhaz"),
               "bhaz(gr = ), has no frmtmb equivalent",
               fixed = TRUE)
  # the pointer to the family is one constant, read by the refusal
  expect_match(frmtmb:::brms_aterm_refusals[["bhaz"]],
               frmtmb:::bhaz_family_pointer, fixed = TRUE)
})
