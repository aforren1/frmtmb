## Punch round 2, NIT 5: what one `$` read on a family costs, against a
## plain list `$` on the same object (unclass()ed), which is the cost a
## family paid before `$.frmtmb_family` existed. The round 0 build is not
## reconstructed: its method is gone from the tree, and the plain list is
## the floor any method is measured from.
##
## Arms, interleaved in one process, each a block of N reads, minimum over
## ROUNDS rounds, Sys.time() (about 2 us resolution; proc.time() ticks at
## 10 ms here):
##   plain    unclass(f)$lpdf, the dispatch-free floor
##   plain2   the same again, the CONTROL; plain2 / plain must read 1.0
##   method   f$lpdf through the installed `$.frmtmb_family`
##   derived  f$link through the installed method
##   bracket  f[["lpdf"]]
##   r1b, r1b_derived  the round 1b method on the same reads (see below)
## Usage: Rscript dev/famlink-p2-dollarcost.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
f <- negbinomial()
p <- unclass(f)
ROUNDS <- 7L
# N is per arm: a plain read is about 100 times cheaper than a method
# read, so one N would leave one arm on clock ticks or the other on
# minutes. Each arm's block is grown until it lasts past 1.2 s.
blk <- function(expr) {
  e <- substitute(expr)
  g <- eval(call("function", as.pairlist(alist(n = )),
                 call("for", as.name("i"), quote(seq_len(n)), e)),
            parent.frame())
  function(n) {
    gc(FALSE)
    t0 <- Sys.time()
    g(n)
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
}
# The method as round 1b shipped it, before the early return for a stored
# name, dispatched on a subclass so it pays the same S3 lookup as the
# installed one. Body copied from R/families.R of that build.
ns <- asNamespace("frmtmb")
`$.famlink_r1b` <- function(x, name) {
  nms <- names(x)
  view <- ns$family_link_view_names(x)
  if (name %in% view) {
    if (name %in% nms) stop("stored link field")
    return(ns$family_link_view(x, name))
  }
  if (name %in% nms) return(.subset2(x, name))
  hit <- c(nms, view)[startsWith(c(nms, view), name)]
  if (length(hit) == 1L) stop("partial match")
  NULL
}
registerS3method("$", "famlink_r1b", `$.famlink_r1b`)
o <- structure(unclass(f), class = c("famlink_r1b", "list"))
stopifnot(identical(o$lpdf, f$lpdf), identical(o$link, f$link))
arms <- list(
  plain = blk(p$lpdf),
  plain2 = blk(p$lpdf),
  method = blk(f$lpdf),
  r1b = blk(o$lpdf),
  derived = blk(f$link),
  r1b_derived = blk(o$link),
  bracket = blk(f[["lpdf"]])
)
N <- sapply(arms, function(a) {
  n <- 1e3
  while (a(n) < 1.2) n <- n * 2
  n
})
cat("reads per block:", paste(names(N), N, collapse = "; "), "\n")
tab <- matrix(NA_real_, ROUNDS, length(arms),
              dimnames = list(NULL, names(arms)))
for (r in seq_len(ROUNDS)) {
  for (a in sample(names(arms))) tab[r, a] <- arms[[a]](N[[a]])
}
us <- 1e6 * apply(tab, 2, min) / N
cat(sprintf("%-8s %8.3f us per read\n", names(us), us), sep = "")
cat(sprintf(paste0("control plain2/plain %.3f; method/plain %.2f; ",
                   "r1b/plain %.2f; r1b/method %.2f; ",
                   "r1b_derived/derived %.2f\n"),
            us[["plain2"]] / us[["plain"]], us[["method"]] / us[["plain"]],
            us[["r1b"]] / us[["plain"]], us[["r1b"]] / us[["method"]],
            us[["r1b_derived"]] / us[["derived"]]))
