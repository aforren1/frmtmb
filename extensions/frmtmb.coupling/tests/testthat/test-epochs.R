## Item 3.2 of dev/extension-gaps-plan.md: a list of epochs of unequal
## length, and a vector pair with a `group` of one label per sample, the
## way frmtmb::frm_periodogram() reads one.
##
## Every call is caught with ep_run() so that on a build without
## `group` the assertions below FAIL rather than stop the block at the
## first error, and each one is seen to fail on the released build.

cls <- "frmtmb_coupling_error"

## A refusal, asserted so that on a build without this code it FAILS
## rather than stopping the block: expect_error(class = ) rethrows an
## error of another class, and every assertion after it would never run.
refuses <- function(expr, regexp) {
  cnd <- tryCatch(expr, error = identity)
  expect_s3_class(cnd, cls)
  expect_match(if (inherits(cnd, "condition")) conditionMessage(cnd) else
                 "", regexp)
}

ep_run <- function(expr) tryCatch(expr, error = identity)

ep_pair <- function(len, alpha = 0.7) {
  x <- lapply(len, stats::rnorm)
  y <- lapply(x, function(v) alpha * v + stats::rnorm(length(v)))
  list(x = x, y = y)
}

## drop the id column so a unit's rows compare with a single-record call
ep_body <- function(d) {
  if (!is.data.frame(d)) return(d)
  d <- d[setdiff(names(d), "id")]
  rownames(d) <- NULL
  d
}

test_that("an epoch list equals the concatenated-and-split call", {
  set.seed(31)
  len <- c(512L, 700L, 333L, 901L)
  p <- ep_pair(len)
  lab <- c("a", "b", "c", "d")
  # one epoch per unit: the list and the vector with a per-sample group
  # describe the same four records, whatever their lengths
  ls <- ep_run(frm_cross_spectrum(stats::setNames(p$x, lab),
                                  stats::setNames(p$y, lab), segments = 4L))
  vs <- ep_run(frm_cross_spectrum(unlist(p$x), unlist(p$y), segments = 4L,
                                  group = rep(lab, len)))
  expect_s3_class(ls, "data.frame")
  expect_identical(ls, vs)
  expect_identical(levels(ls$id), lab)
  # and each unit is exactly the call on that epoch alone
  for (i in seq_along(lab)) {
    one <- frm_cross_spectrum(p$x[[i]], p$y[[i]], segments = 4L)
    got <- if (is.data.frame(ls)) ls[ls$id %in% lab[i], , drop = FALSE]
    expect_identical(ep_body(got), one, label = lab[i])
  }
})

test_that("pooled epochs equal the split call when the lengths agree", {
  set.seed(32)
  # two units of two epochs each, every epoch 512 samples; segments = 4
  # makes the segment 256 long, which divides 512, so the concatenated
  # record and the pooled epochs cut the same blocks
  p <- ep_pair(rep(512L, 4L))
  g <- c("s1", "s1", "s2", "s2")
  ls <- ep_run(frm_cross_spectrum(p$x, p$y, segments = 4L, group = g))
  vs <- ep_run(frm_cross_spectrum(unlist(p$x), unlist(p$y), segments = 4L,
                                  group = rep(g, each = 512L)))
  expect_s3_class(ls, "data.frame")
  expect_identical(ls, vs)
  expect_identical(unique(ls$n), 4L)
})

test_that("pooled epochs of unequal length are clean spans of one record", {
  set.seed(33)
  # a unit's epochs are read as the clean spans of one record, so the
  # list must equal the single-record call with an NA between epochs,
  # which is the gap machinery this package already ships
  len <- c(700L, 1043L, 699L)
  p <- ep_pair(len)
  gap_x <- unlist(lapply(p$x, function(v) c(v, NA)))
  gap_y <- unlist(lapply(p$y, function(v) c(v, NA)))
  for (seg in c(4L, 8L, 16L)) {
    ls <- ep_run(frm_cross_spectrum(p$x, p$y, segments = seg,
                                    group = rep("u", 3L)))
    ref <- frm_cross_spectrum(gap_x, gap_y, segments = seg)
    expect_identical(ep_body(ls), ref, label = paste("segments", seg))
  }
  # the concatenation WITHOUT the gaps is a different record: segments
  # there cross the epoch boundaries, so the two must disagree
  joined <- frm_cross_spectrum(unlist(p$x), unlist(p$y), segments = 8L)
  pooled <- ep_run(frm_cross_spectrum(p$x, p$y, segments = 8L,
                                      group = rep("u", 3L)))
  expect_false(isTRUE(all.equal(ep_body(pooled), joined)))
  # and no epoch gives a segment it cannot hold: 2441 samples at 8
  # segments is 305 a segment, and the epochs hold 2, 3 and 2
  expect_identical(unique(pooled$n), 7L)
})

test_that("group labels name the units in order of first appearance", {
  set.seed(34)
  p <- ep_pair(c(400L, 420L, 440L))
  xs <- ep_run(frm_cross_spectrum(p$x, p$y, segments = 4L,
                                  group = c("z", "a", "z")))
  expect_identical(levels(xs$id), c("z", "a"))
  # unnamed epochs with no group are numbered
  xs2 <- ep_run(frm_cross_spectrum(p$x, p$y, segments = 2L))
  expect_identical(levels(xs2$id), c("1", "2", "3"))
})

test_that("the list and group forms refuse what they cannot read", {
  set.seed(35)
  p <- ep_pair(c(300L, 400L))
  refuses(frm_cross_spectrum(p$x, p$y, group = "a"),
               "one label per epoch")
  refuses(frm_cross_spectrum(p$x, p$y, group = c("a", NA)),
               "has NA")
  refuses(frm_cross_spectrum(p$x, p$y[1L]),
               "same number of epochs")
  bad <- p$y
  bad[[2L]] <- bad[[2L]][-1L]
  refuses(frm_cross_spectrum(p$x, bad),
               "epoch 2 has 400 samples in `x` and 399")
  refuses(frm_cross_spectrum(stats::setNames(p$x, c("a", "b")),
                                  stats::setNames(p$y, c("b", "a"))),
               "name their epochs differently")
  refuses(frm_cross_spectrum(p$x, p$y[[1L]]),
               "both be lists of epochs")
  refuses(frm_cross_spectrum(list(p$x), list(p$y)),
               "is not a vector")
  refuses(frm_cross_spectrum(list(letters), list(letters)),
               "x\\[\\[1\\]\\]")
  v <- stats::rnorm(600)
  refuses(frm_cross_spectrum(v, v, group = rep("a", 599)),
               "one label per sample")
  refuses(frm_cross_spectrum(v, v, group = c(rep("a", 599), NA)),
               "has NA")
  m <- matrix(stats::rnorm(1200), 600)
  refuses(frm_cross_spectrum(m, m, group = c("a", "b")),
               "already the units")
  refuses(frm_cross_spectrum(data.frame(a = v), data.frame(a = v)),
               "data frame")
})
