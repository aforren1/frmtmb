# What the new export actually exposes, what the Rd page says, and
# whether the two agree.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))

cat("== 1. what is exported ==\n")
ex <- getNamespaceExports("frmtmb")
for (n in c("frm_install_generics", "frm_bind_generic",
            "frm_adopt_target", "frm_is_generic", "frm_generic_owners"))
  cat(sprintf("  %-22s exported: %-5s exists in namespace: %s\n", n,
              n %in% ex, exists(n, envir = asNamespace("frmtmb"),
                                inherits = FALSE)))
cat("  frmtmb.sample re-exports it: ",
    "frm_install_generics" %in%
      getNamespaceExports("frmtmb.sample"), "\n")
cat("  reachable as frmtmb::frm_install_generics: ",
    is.function(frmtmb::frm_install_generics), "\n")
cat("  formals: ", paste(deparse(args(frmtmb::frm_install_generics)),
                         collapse = " "), "\n\n")

cat("== 2. the Rd page, RENDERED, not read from source ==\n")
rd <- tools::Rd_db("frmtmb")[["frmtmb-sampling-api.Rd"]]
f <- tempfile()
tools::Rd2txt(rd, out = f, options = list(underline_titles = FALSE))
txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
want <- c("The borrowed-generic seam",
          "frm_install_generics(pkgname, owners)",
          "S3method(pkg::generic, class)",
          "bare `UseMethod()`",
          "Call it FIRST in the extension's own",
          "frm_install_generics")
for (w in want)
  cat(sprintf("  %-46s present: %s\n", substr(w, 1, 44),
              grepl(w, txt, fixed = TRUE)))
cat("  the page lists frm_bind_generic (it should NOT): ",
    grepl("frm_bind_generic", txt, fixed = TRUE), "\n")
cat("  stray percent signs in the section (Rd comment hazard): ",
    length(grep("%", strsplit(txt, "\n")[[1]], value = TRUE)), "\n")
cat("  usage/alias block has frm_install_generics: ",
    "frm_install_generics" %in% unlist(lapply(rd, function(x)
      if (identical(attr(x, "Rd_tag"), "\\alias"))
        trimws(paste(unlist(x), collapse = "")))), "\n")
sec <- regmatches(txt, regexpr(
  "The borrowed-generic seam(.|\n)*?\n\n\n", txt))
cat("\n--- the rendered section ---\n")
cat(substr(sec, 1, 2200), "\n")

cat("\n== 3. the new export called on a SEALED namespace ==\n")
r <- tryCatch(frmtmb::frm_install_generics("stats",
                 owners = list(sd = "stats")),
              error = function(e) conditionMessage(e))
cat("  frm_install_generics(\"stats\", list(sd = \"stats\")): ", r, "\n")
r2 <- tryCatch(frmtmb::frm_install_generics("nosuchpkg",
                 owners = list(a = "b")),
               error = function(e) conditionMessage(e))
cat("  a package that is not loaded: ", r2, "\n")
r3 <- tryCatch(frmtmb::frm_install_generics("frmtmb",
                 owners = list(zzz_not_a_generic = "stats")),
               error = function(e) conditionMessage(e))
cat("  a name frmtmb does not have, owner without it: ",
    paste(r3, collapse = " "), "\n")
cat("  after that call, is zzz_not_a_generic bound in frmtmb? ",
    exists("zzz_not_a_generic", envir = asNamespace("frmtmb"),
           inherits = FALSE), "\n")

cat("\n== 4. bbmle's non-generic parnames ==\n")
q(library(frmtmb.sample))
q(library(bbmle))
cat("  search: ", paste(head(search(), 4), collapse = " "), "\n")
g <- get("parnames", envir = globalenv())
cat("  parnames from globalenv resolves to: ",
    environmentName(environment(g)), "\n")
cat("  bbmle::parnames still reachable with ::: ",
    is.function(bbmle::parnames), " formals ",
    paste(deparse(args(bbmle::parnames)), collapse = " "), "\n")
o <- structure(numeric(3), names = NULL)
cat("  bbmle::parnames(o) works: ",
    tryCatch({ bbmle::parnames(o); "yes" },
             error = function(e) conditionMessage(e)), "\n")
cat("  plain parnames(o) (the masked call): ",
    tryCatch({ print(parnames(o)); "returned" },
             error = function(e) substr(conditionMessage(e), 1, 70)), "\n")
cat("  and bbmle's REPLACEMENT form, which :: cannot spell: ",
    exists("parnames<-", envir = asNamespace("bbmle")), "\n")
r4 <- tryCatch({ x <- numeric(3); `parnames<-`(x, c("a","b","c"))
  "assigned" }, error = function(e) substr(conditionMessage(e), 1, 70))
cat("  `parnames<-`(x, ...) from the global environment: ", r4, "\n")
cat("DONE\n")
