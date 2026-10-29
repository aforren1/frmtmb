# Review of lane meme, short check of fix round 1: shapes for the
# sampling-route table comparison. Sourced after review-meme-shapes.R.
sc_data <- function() {
  d <- rv_data()
  set.seed(31)
  d$th <- factor(rep(c("a", "b"), length.out = nrow(d)))
  d$yo <- sample(1:4, nrow(d), TRUE)
  d$u <- rnorm(nrow(d))
  d
}
sc_shapes <- function(bf, mvbf, set_rescor, cumulative, gp) {
  list(
    mixed3 = list(f = bf(y ~ me(xg, sg, gr = g) + me(z, sz) +
                           me(wg, sg, gr = g))),
    mixed4 = list(f = bf(y ~ me(xg, sg, gr = g) + me(z, sz) +
                           me(wg, sg, gr = g) + me(x, sx))),
    reme = list(f = bf(y ~ x + (1 + x | g) + me(xg, sg, gr = g) +
                         me(wg, sg, gr = g))),
    mvre = list(f = mvbf(bf(y ~ x + (1 | g)), bf(y2 ~ x + (1 + x | g))) +
                  set_rescor(TRUE)),
    dparre = list(f = bf(y ~ x + (1 | g), sigma ~ (1 | h))),
    thres = list(f = bf(yo | thres(gr = th) ~ x), fam = cumulative()),
    thres_eq = list(f = bf(yo | thres(gr = th) ~ x),
                    fam = cumulative(threshold = "equidistant")),
    nosdx = list(f = bf(y ~ me(x) + z)),
    gp = list(f = bf(y ~ gp(u) + x)))
}
