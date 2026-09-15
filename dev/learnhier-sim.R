# Lane `learnhier` (item 2.2): the two designs, and the truth each one
# has to be scored against.
#
# WHY THE TRUTH IS NOT ONE CONSTANT VECTOR. Item 2.2 asks whether a
# correlated block over EVERY parameter recovers its correlations. A
# design that draws independent deviations cannot answer that: it only
# says whether the estimator manufactures a correlation that is not
# there. So both designs below draw a block with real off-diagonal
# structure, and different signs, so that shrinkage toward zero and a
# sign error are separately visible.
#
# For `bandit2arm_delta()` the block the model fits IS the block the
# design drew: `alpha` is estimated on the logit and `tau` on the log,
# which are the scales the deviations are drawn on, so the drawn
# correlation is the truth and coverage against it is a coverage.
#
# For `rlddm()` under `ndt_group(id)` it is NOT. The model estimates
# `qlogis(ndt_i / floor_i)`, where `floor_i` is that learner's own
# fastest response, which is an observed minimum of the learner's own
# drawn data. A learner with a wider boundary is slower, so its floor is
# later, so it needs a SMALLER fraction to express the same non-decision
# time: the `bs` and `ndt` deviations are correlated in the fitted
# parameterization before any fit. lh_rlddm_data() therefore records the
# INDUCED truth, computed from the drawn truths and the realized floors
# with no fit at all, alongside the drawn one.
#
# A floor is an observed minimum, so it moves with the trial count.
# Every figure this lane reports says which count it came from; the
# designs here are the plan's realistic scale, 100 learners by 200
# trials.

# ------------------------------------------------------------- truths

# Correlations chosen so that no two pairs share a value and the signs
# differ: an estimator that returned the average of the block, or that
# flipped a sign, is visible in the table rather than only in a norm.
lh_bandit_truth <- list(
  alpha = 0.35, tau = 3,
  sd = c(alpha = 0.5, tau = 0.3),
  R = matrix(c(1, 0.5,
               0.5, 1), 2L, 2L,
             dimnames = list(c("alpha", "tau"), c("alpha", "tau"))))

lh_rlddm_truth <- list(
  alpha = 0.35, drift = 2.5, bs = 1.5, ndt = 0.25, bias = 0.5,
  sd = c(alpha = 0.5, drift = 1.0, bs = 0.2, ndt = 0.15),
  R = matrix(c(1.0,  0.4,  0.0,  0.0,
               0.4,  1.0, -0.3,  0.0,
               0.0, -0.3,  1.0,  0.3,
               0.0,  0.0,  0.3,  1.0), 4L, 4L, byrow = TRUE,
             dimnames = list(c("alpha", "drift", "bs", "ndt"),
                             c("alpha", "drift", "bs", "ndt"))))

# The zero-correlation arm of the bandit design, kept because it answers
# the other half of the question: an estimator that recovers a real
# correlation is not much use if it also reports one that is not there.
lh_bandit_truth0 <- within(lh_bandit_truth, R <- diag(2))

# ------------------------------------------------------- the draw

# One draw of `n` correlated deviations with the stated standard
# deviations. chol() is used rather than eigen() so that a truth matrix
# that is not positive definite is an error here rather than a silent
# reflection of a negative eigenvalue.
lh_draw_dev <- function(n, sd, R) {
  L <- chol(R)
  Z <- matrix(stats::rnorm(n * ncol(R)), n, ncol(R))
  D <- Z %*% L
  D <- sweep(D, 2L, sd, "*")
  colnames(D) <- names(sd)
  D
}

# 100 learners x 200 trials of a stationary two-armed bandit, with a
# correlated deviation on the learning rate's logit and on the inverse
# temperature's log.
lh_bandit_data <- function(seed, ns = 100L, nt = 200L,
                           truth = lh_bandit_truth) {
  d <- frmtmb.learn::frm_task_design("bandit2arm", n_subject = ns,
                                     n_trial = nt, seed = seed)
  set.seed(seed + 1L)
  U <- lh_draw_dev(ns, truth$sd, truth$R)
  i <- as.integer(d$id)
  d$choice <- frmtmb.learn::frm_task_simulate(
    frmtmb.learn::bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(stats::qlogis(truth$alpha) +
                                        U[i, "alpha"]),
                tau = exp(log(truth$tau) + U[i, "tau"])),
    seed = seed)[[1L]]$choice
  attr(d, "dev_drawn") <- U
  # for this family the fitted block IS the drawn block
  attr(d, "dev_fitted") <- U
  d
}

# The same task frame under a diffusion choice rule. Four correlated
# deviations, `bias` held at 0.5.
lh_rlddm_data <- function(seed, ns = 100L, nt = 200L,
                          truth = lh_rlddm_truth) {
  d <- frmtmb.learn::frm_task_design("bandit2arm", n_subject = ns,
                                     n_trial = nt, seed = seed)
  set.seed(seed + 2L)
  U <- lh_draw_dev(ns, truth$sd, truth$R)
  i <- as.integer(d$id)
  ndt_i <- truth$ndt * exp(U[, "ndt"])
  s <- frmtmb.learn::frm_task_simulate(
    frmtmb.learn::rlddm(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(stats::qlogis(truth$alpha) +
                                        U[i, "alpha"]),
                drift = truth$drift + U[i, "drift"],
                bs = truth$bs * exp(U[i, "bs"]),
                ndt = ndt_i[i], bias = truth$bias),
    seed = seed)[[1L]]
  fl <- as.numeric(tapply(s$rt, s$id, min))[match(levels(s$id),
                                                  levels(s$id))]
  attr(s, "ndt_subject") <- ndt_i
  attr(s, "own_floor") <- fl
  attr(s, "dev_drawn") <- U
  # THE BLOCK THE MODEL ACTUALLY FITS. Three columns are the drawn
  # deviations; the fourth is the logit of each learner's realized
  # fraction of its own floor, which is what `ndt ~ (1 | p | id)`
  # estimates under ndt_group(). No fit enters this.
  attr(s, "dev_fitted") <- cbind(
    alpha = U[, "alpha"], drift = U[, "drift"], bs = U[, "bs"],
    ndt = stats::qlogis(pmin(pmax(ndt_i / fl, 1e-8), 1 - 1e-8)))
  s
}

# ------------------------------------------- reading a block off a fit

# The lower-triangle pair labels of a d x d block, in the order
# confint_varcorr() reports them, so a table row can be matched by name
# rather than by position.
lh_pairs <- function(nms) {
  d <- length(nms)
  ij <- which(lower.tri(diag(d)), arr.ind = TRUE)
  ij <- ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
  list(lab = paste0(nms[ij[, "col"]], "~", nms[ij[, "row"]]),
       i = ij[, "col"], j = ij[, "row"])
}

# The realized correlations and standard deviations of a deviation
# matrix, named the way lh_pairs() names them.
lh_block_stats <- function(D) {
  nms <- colnames(D)
  pr <- lh_pairs(nms)
  cr <- stats::cor(D)
  c(stats::setNames(apply(D, 2L, stats::sd), paste0("sd_", nms)),
    stats::setNames(cr[cbind(pr$j, pr$i)], paste0("cor_", pr$lab)))
}
