# Build the frmtmb documentation site and the four extension subsites.
#
#   Rscript dev/build-docs.R                 # all five, core first
#   Rscript dev/build-docs.R frmtmb.sample   # one package
#
# Run it from the repository root. The core site lands in docs/ and
# each extension lands in docs/<package>/, so one GitHub Pages
# deployment serves all five. There is no site workflow in
# .github/workflows: docs/ is checked in, and the site is built here
# and committed, which is the practice this repository already
# follows.
#
# WHY THE INSTALL COMES FIRST. The point of one site is that
# `frmtmb.sample::frm_sample()` in the core prose, and
# `frmtmb::get_prior()` in the extension prose, come out as links.
# downlit does that resolution, and it needs two things from the
# package being referenced: the package installed, so that it can look
# the topic up in the help index and learn which .Rd file holds it;
# and the address of that package's site. So every package is
# installed before any site is built, and a stale install of one
# package silently degrades the links in the other four.
#
# WHY THE SITE METADATA IS WRITTEN BY HAND. downlit reads the site
# address from `system.file("pkgdown.yml", package = )` and, failing
# that, fetches <DESCRIPTION URL>/pkgdown.yml over the network. pkgdown
# writes that file into the built site, not into the installed package,
# so neither source is available for a subsite that has never been
# deployed. write_site_metadata() below puts the file where downlit
# looks for it. It is derived from the same _pkgdown.yml and
# vignettes/ that the build reads, and check_metadata() compares it
# against what pkgdown wrote, so the two cannot drift apart unnoticed.
#
# COST. A full run takes about eleven minutes on the machine this was
# written on, of which the core site is about eight: every vignette in
# every package is rebuilt, and the core vignettes fit models. The five
# installs are about half a minute of that. Name one package on the
# command line to rebuild only its site; the installs still run,
# because the links in that one site point at the other four.

PKGS <- c(
  frmtmb          = ".",
  frmtmb.sample   = "extensions/frmtmb.sample",
  frmtmb.latent   = "extensions/frmtmb.latent",
  frmtmb.ode      = "extensions/frmtmb.ode",
  frmtmb.ddm      = "extensions/frmtmb.ddm"
)

if (!file.exists("DESCRIPTION") ||
      read.dcf("DESCRIPTION", "Package")[[1]] != "frmtmb") {
  stop("run this from the repository root")
}

# Installing into .libPaths()[1] replaces whatever build of these five
# packages is already there. Set FRMTMB_DOCS_LIB to keep the working
# library untouched.
LIB <- Sys.getenv("FRMTMB_DOCS_LIB", unset = .libPaths()[1])
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))

# pkgdown renders vignettes in a child process, which reads the
# library path from the environment rather than from this session.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

`%||%` <- function(x, y) if (is.null(x)) y else x

site_url <- function(dir) {
  yaml::read_yaml(file.path(dir, "_pkgdown.yml"))$url
}

vignette_names <- function(dir) {
  files <- list.files(file.path(dir, "vignettes"), pattern = "[.]Rmd$")
  files <- files[!startsWith(files, "_")]
  stats::setNames(sub("[.]Rmd$", ".html", files), sub("[.]Rmd$", "", files))
}

write_site_metadata <- function(pkg, dir) {
  url <- site_url(dir)
  meta <- list(
    articles = as.list(vignette_names(dir)),
    urls = list(
      reference = paste0(url, "/reference"),
      article = paste0(url, "/articles")
    )
  )
  yaml::write_yaml(meta, file.path(LIB, pkg, "pkgdown.yml"))
}

# The primed metadata is a prediction of what pkgdown will write. If
# the prediction is wrong the links point at pages that do not exist,
# and nothing else in the build would say so.
check_metadata <- function(pkg, dir) {
  dest <- yaml::read_yaml(file.path(dir, "_pkgdown.yml"))$destination %||% "docs"
  built <- file.path(dir, dest, "pkgdown.yml")
  if (!file.exists(built)) return(invisible())
  a <- yaml::read_yaml(file.path(LIB, pkg, "pkgdown.yml"))
  b <- yaml::read_yaml(built)
  if (!identical(a$urls, b$urls) || !identical(a$articles, b$articles)) {
    warning(pkg, ": primed site metadata disagrees with the built site",
            call. = FALSE, immediate. = TRUE)
  }
  invisible()
}

wanted <- commandArgs(trailingOnly = TRUE)
if (length(wanted) == 0) wanted <- names(PKGS)
if (!all(wanted %in% names(PKGS))) {
  stop("unknown package: ", paste(setdiff(wanted, names(PKGS)), collapse = ", "))
}

started <- Sys.time()

# Help pages are what downlit reads, so --no-docs would break the
# links this whole script exists to produce. Vignettes cost nothing
# here: R CMD INSTALL of a source directory does not build them, and
# pkgdown renders its own copies from vignettes/ anyway.
for (pkg in names(PKGS)) {
  message("--- installing ", pkg)
  utils::install.packages(
    PKGS[[pkg]], lib = LIB, repos = NULL, type = "source",
    INSTALL_opts = "--no-multiarch"
  )
  write_site_metadata(pkg, PKGS[[pkg]])
}

installed <- Sys.time()
message(sprintf("--- installed five packages in %.1f min",
                as.numeric(difftime(installed, started, units = "mins"))))

for (pkg in wanted) {
  message("--- building ", pkg)
  pkgdown::build_site(
    pkg = PKGS[[pkg]],
    preview = FALSE,
    install = FALSE,
    new_process = FALSE
  )
  check_metadata(pkg, PKGS[[pkg]])
}

message(sprintf("--- done in %.1f min",
                as.numeric(difftime(Sys.time(), started, units = "mins"))))
