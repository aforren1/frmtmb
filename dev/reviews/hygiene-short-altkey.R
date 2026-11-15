# Short check of lane hygiene's fix round, F1: a candidate key that
# never collapses the cases hygiene-short-regkey.R collapses and still
# takes a reload for the same registration. Not the lane's code.
#   Rscript dev/reviews/hygiene-short-altkey.R
# Rule: walk from the function's environment up to the first named one
# (namespace, global, package); serialize every frame's bindings; any
# reference object met on the way (an unnamed environment, an external
# pointer: the refhook sees both) makes the key unique.
alt_key <- function(fn) {
  env <- environment(fn) %||% baseenv()
  named <- function(e) isNamespace(e) || identical(e, globalenv()) ||
    identical(e, baseenv()) || identical(e, emptyenv()) ||
    nzchar(environmentName(e))
  frames <- list()
  ref <- FALSE
  while (!named(env)) {
    frames[[length(frames) + 1L]] <- tryCatch(
      serialize(as.list(env, all.names = TRUE, sorted = TRUE), NULL,
                refhook = function(x) { ref <<- TRUE; NULL }),
      error = function(e) { ref <<- TRUE; NULL })
    env <- parent.env(env)
  }
  if (ref) return(new.env())
  list(owner = environmentName(topenv(environment(fn) %||% baseenv())),
       code = paste(deparse(fn), collapse = "\n"), frames = frames)
}
same <- function(f, g) identical(alt_key(f), alt_key(g))

outer <- function(mode) { force(mode); function(label) { force(label)
  function(spec, frame) if (mode == "refuse") stop("x") } }
via_local <- function(mode) { force(mode)
  local(function(spec, frame) if (mode == "refuse") stop("x")) }
mk_cfg <- function(cfg) { force(cfg); function(spec, frame) cfg$refuse }
mk_ptr <- function(ptr) { force(ptr); function(spec, frame) ptr }
p1 <- get("C_Cdqrls", asNamespace("stats"))$address
p2 <- get("C_DoubleCentre", asNamespace("stats"))$address
make_check <- function(tag) { force(tag); function(spec, frame) tag }
# sample's provider: a closure made in .onLoad(libname, pkgname)
onload <- function(libname, pkgname) function(spec, frame) pkgname
environment(onload) <- asNamespace("stats")
cat("ALT nested, different outer value, same:",
    same(outer("refuse")("l"), outer("log")("l")), "(want FALSE)\n")
cat("ALT local() in a factory, same:",
    same(via_local("refuse"), via_local("log")), "(want FALSE)\n")
cat("ALT two empty environments, same:",
    same(mk_cfg(new.env()), mk_cfg(new.env())), "(want FALSE)\n")
cat("ALT two external pointers, same:",
    same(mk_ptr(p1), mk_ptr(p2)), "(want FALSE)\n")
cat("ALT factory, same value twice, same:",
    same(make_check("a"), make_check("a")), "(want TRUE)\n")
cat("ALT factory, different values, same:",
    same(make_check("a"), make_check("b")), "(want FALSE)\n")
cat("ALT nested, same values twice, same:",
    same(outer("log")("l"), outer("log")("l")), "(want TRUE)\n")
cat("ALT .onLoad closure, same call twice, same:",
    same(onload("L", "p"), onload("L", "p")), "(want TRUE)\n")
cat("ALT namespace function twice, same:",
    same(stats::model.frame, stats::model.frame), "(want TRUE)\n")
