# Reinforcement learning: a sequential family from scratch

frmtmb has no reinforcement-learning family. This page writes one and
fits it hierarchically through the ordinary formula grammar. The whole
family is under 200 lines, the taped likelihood is 23 of them, and about
a third of the rest is the text of the messages it owes a user when it
refuses something.

The point is not the model. The point is the seam. A likelihood that
reads a whole sequence rather than one row at a time reaches the core
through
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md),
and everything else keeps working: random effects, links, standard
errors, prediction.

Read
[`vignette("frmtmb")`](https://aforren1.github.io/frmtmb/articles/frmtmb.md)
first for the grammar, and
[`?frmtmb_structure`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the reference description of the protocol used here.

## The model

The task is a two-armed bandit. On each trial a subject picks one of two
arms and receives a payoff. The Rescorla-Wagner delta rule says the
subject keeps a value estimate `Q` for each arm, starts both at zero,
and after each choice moves the chosen arm’s estimate toward what it
just paid:

    Q[chosen] <- Q[chosen] + alpha * (reward - Q[chosen])

The choice itself is a softmax over the two values, which for two arms
is a logistic function of their difference:

    P(arm 1) <- plogis(beta * (Q1 - Q2))

`alpha` is the learning rate, on `(0, 1)`. `beta` is the inverse
temperature, on `(0, Inf)`: a large `beta` means the subject almost
always takes the arm it currently values more. This is the model
hBayesDM calls `bandit2arm_delta` (Ahn, Haines and Zhang 2017,
*Computational Psychiatry* 1, 24). None of hBayesDM’s code is used here;
the equations above are the published model, and the family below is
written from them.

Each of the two parameters becomes a distributional parameter with its
own linear predictor, so a learning rate can carry a condition effect
and a subject-level random intercept exactly the way a mean would.

## Why this needs the structured protocol

The likelihood does factorize over trials: given the parameters, `Q` is
a deterministic function of the subject’s own earlier choices and
payoffs, so each trial contributes one Bernoulli factor.

What it is not is ROWWISE. A row’s factor cannot be computed from that
row’s parameters and response alone, because `Q` carries the whole
history. `frmtmb_family(lpdf =)` is handed one row’s worth of
everything, so it cannot express this. `frmtmb_structure(loglik =)` is
handed the whole response and returns one number, so it can.

## The family

### The data a trial carries

``` r

# reward() carries what each arm WOULD have paid on this trial, in arm
# order. The likelihood only ever reads the chosen arm's entry, so data
# that records the received outcome alone passes it twice. The second
# column is what makes the simulator coherent, because a simulated
# choice needs the payoff of the arm the subject did not take in the
# data. vreal(pay1, pay2) would carry the same two columns. The registry
# is what lets a family name its data in the words of its own
# literature, and `reward` is the word this literature uses.
frmtmb_register_aterm("reward", arity = 2L)
```

### The block: from a long data frame to subject by trial

The recursion is sequential in trials and independent across subjects,
so it should walk trials once and update every subject at each step.
That needs the rows arranged as a subject-by-trial matrix, which is what
the `frame_block` slot is for. It runs once, at frame assembly, and
everything it returns is plain data that gets saved inside the fit.

``` r

# The frame block is DATA, resolved once, and it is what turns a long
# data frame into the subject-by-trial layout the recursion walks. Every
# subject gets a row of `idx` holding its row numbers in trial order,
# padded on the right to the longest subject with a repeat of its own
# first row. `mask` is 1 on a real trial and 0 on a pad. It multiplies
# both the log-likelihood and the learning rate, so a padded cell
# contributes nothing and updates nothing.
rw_block <- function(resp, spec, av, mf, y, n) {
  rw <- resp$family[["rw"]]
  gv <- eval(rw[["subject_expr"]], mf, resp$formula_env)
  if (anyNA(gv)) {
    stop("rw_delta(): every row needs a subject, and ",
         deparse1(rw[["subject_expr"]]), " has ", sum(is.na(gv)),
         " missing value(s)", call. = FALSE)
  }
  gv <- factor(gv)
  tv <- if (is.null(rw[["trial_expr"]])) {
    # no trial column: the data frame's own row order within a subject
    out <- integer(n)
    for (g in split(seq_len(n), gv)) out[g] <- seq_along(g)
    out
  } else {
    v <- eval(rw[["trial_expr"]], mf, resp$formula_env)
    if (anyNA(v)) {
      stop("rw_delta(): the trial variable '", deparse1(rw[["trial_expr"]]),
           "' has missing values, so the order of the recursion is ",
           "undefined at those rows", call. = FALSE)
    }
    as.numeric(v)
  }
  rows <- lapply(split(seq_len(n), gv), function(r) r[order(tv[r])])
  key <- paste(as.integer(gv), tv, sep = "|")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1L]
    stop("rw_delta(): trial numbers must be unique within a subject; ",
         sum(key == dup), " rows share one. A delta rule updates once ",
         "per trial, so two rows at one trial have no order to learn in",
         call. = FALSE)
  }
  len <- lengths(rows)
  nt <- max(len)
  # The matrix() before each transpose is load-bearing: vapply() drops to
  # a plain vector when nt == 1, and t() of a vector is 1-by-n_subj, so
  # the block would come back with every subject read as one more trial
  # of a single learner and Q carried across subject boundaries.
  idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], nt - length(r))),
                         integer(nt)), nrow = nt))
  mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(nt) <= l),
                          numeric(nt)), nrow = nt))
  list(idx = idx, mask = mask, len = len, n_subj = length(rows),
       n_trial = nt, subject = gv, trial = tv, levels = levels(gv), n = n)
}
```

Two things to notice. `idx` holds row NUMBERS, not values, so the
recursion can index any per-row quantity with it: the response, the
payoffs, and the linear predictors of `alpha` and `beta`. And subjects
with different trial counts are padded to the longest one, with `mask`
marking the padding. The mask is what makes the padding free of branches
on the tape: it multiplies the log-likelihood, so a padded cell
contributes nothing, and it multiplies the learning rate, so a padded
cell learns nothing.

Padding is not the only way.
[`frmtmb.latent::hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.html)
sorts its sequences by decreasing length instead, which makes the
sequences still running at step `s` a prefix of the order and needs no
padded cells at all. That is faster and harder to read. Measured on 40
subjects whose trial counts run from 30 to 100, the padded design
carries 41 percent empty cells and builds its tape in 0.09 seconds
against the balanced design’s 0.08, so for a tutorial the mask wins.

### The likelihood

``` r

# Broadcast a dpar the core handed back as a single value. rep(v, each =)
# strips the advector class (see the RTMB notes in dev/); multiplying by
# a vector of ones does not.
rw_bcast <- function(v, n) if (length(v) == 1L) v * rep(1, n) else v

# The taped log-likelihood of the whole response, as one AD scalar.
#
# One iteration per TRIAL, each one a handful of vector operations over
# all subjects at once. Writing it the other way round, one iteration
# per row with Q[s] <- ... inside, is the shape RTMB punishes: an
# element-wise assignment into a taped vector costs about three orders
# of magnitude more than the vector operation it replaces.
rw_loglik <- function(y, dpars, aterms, weights, block, extra) {
  # weights() is refused in check_spec, so `weights` is 1 here. A trial's
  # factor cannot be reweighted on its own: its value depends on every
  # earlier trial of the same subject.
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  n <- block[["n"]]
  alpha <- rw_bcast(dpars[["alpha"]], n)
  beta <- rw_bcast(dpars[["beta"]], n)
  pay1 <- aterms[["reward1"]]
  pay2 <- aterms[["reward2"]]
  q1 <- rep(0, nrow(idx))
  q2 <- q1
  ll <- 0
  for (t in seq_len(ncol(idx))) {
    i <- idx[, t]
    m <- msk[, t]
    c1 <- y[i]
    eta <- beta[i] * (q1 - q2)
    # log P(choice) = c1 * eta - log(1 + exp(eta)), folded with
    # logspace_add so a decisive subject does not overflow
    ll <- ll + sum(m * (c1 * eta - RTMB::logspace_add(0 * eta, eta)))
    # the mask enters through the learning rate, so a padded cell learns
    # nothing and no branch reaches the tape
    a <- alpha[i] * m
    q1 <- q1 + (a * c1) * (pay1[i] - q1)
    q2 <- q2 + (a * (1 - c1)) * (pay2[i] - q2)
  }
  ll
}
```

The loop body is the delta rule, written for every subject at once. `q1`
and `q2` are vectors with one entry per subject, and one iteration of
the loop advances all of them. There is no conditional anywhere: the
chosen arm is selected by multiplying by the response, which is data,
and the padding is selected by multiplying by the mask, which is also
data.

`weights` arrives and is ignored, which the family is allowed to do only
because it refuses [`weights()`](https://rdrr.io/r/stats/weights.html)
outright further down. A row weight has no meaning here. A trial’s
factor depends on every earlier trial of the same subject, so there is
nothing separable to reweight.

### Fitted values and draws

The same recursion in plain numbers answers two more questions.

``` r

# The same recursion off the tape, in plain numbers, at the estimates.
# It answers two questions: the per-trial choice probability, which is
# the fitted value, and a forward draw, which is the simulator. They
# share one function because a draw is a replay whose choices come from
# rbinom() instead of from the data.
rw_replay <- function(dpars, aterms, block, y, draw = FALSE) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  n <- block[["n"]]
  alpha <- rep(as.numeric(dpars[["alpha"]]), length.out = n)
  beta <- rep(as.numeric(dpars[["beta"]]), length.out = n)
  pay1 <- aterms[["reward1"]]
  pay2 <- aterms[["reward2"]]
  p <- rep(NA_real_, n)
  q1 <- rep(0, nrow(idx))
  q2 <- q1
  for (t in seq_len(ncol(idx))) {
    i <- idx[, t]
    m <- msk[, t]
    keep <- m == 1
    pt <- stats::plogis(beta[i] * (q1 - q2))
    p[i[keep]] <- pt[keep]
    c1 <- if (draw) stats::rbinom(length(i), 1L, pt) else y[i]
    if (draw) y[i[keep]] <- c1[keep]
    a <- alpha[i] * m
    q1 <- q1 + (a * c1) * (pay1[i] - q1)
    q2 <- q2 + (a * (1 - c1)) * (pay2[i] - q2)
  }
  list(p = p, y = y)
}

# The three post-fit slots, each one line of work on top of the replay.
rw_parts <- function(fit) {
  rspec <- single_response(fit, "a rw_delta() fit")
  list(dpars = eval_dpars(fit)[[rspec$resp_name]],
       aterms = fit$frame[["aterm_values"]][[rspec$resp_name]],
       y = fit$frame[["y"]][[rspec$resp_name]])
}

rw_fitted <- function(fit, block) {
  pp <- rw_parts(fit)
  rw_replay(pp$dpars, pp$aterms, block, pp$y)$p
}
```

The per-trial choice probability is the natural fitted value for this
family, and its binomial variance is what pearson residuals divide by. A
forward draw is the same walk with
[`rbinom()`](https://rdrr.io/r/stats/Binomial.html) in place of the
observed choice, which is why one function serves both.

``` r

# The structured simulator. simulate(), posterior_predict() and
# frm_simulate() all reach it through the same context, so the draw is
# written once. It walks each subject forward, drawing a choice from the
# current Q pair and learning from the payoff that choice earns.
rw_sim <- function(ctx) {
  blk <- ctx[["block"]]
  rw_replay(ctx[["dpars"]], ctx[["aterms"]], blk,
            y = rep(0, blk[["n"]]), draw = TRUE)$y
}
```

### The family object

``` r

# The family. Two distributional parameters, each with its own linear
# predictor: alpha on (0, 1) through a logit link, beta on (0, Inf)
# through a log link. `primary_dpars = "alpha"` sends the main
# right-hand side of the formula to the learning rate.
rw_delta <- function(subject, trial = NULL) {
  subject_expr <- substitute(subject)
  trial_expr <- substitute(trial)
  if (is.null(subject_expr)) {
    stop("rw_delta(subject =) names the column that separates one ",
         "learner's trial sequence from the next", call. = FALSE)
  }
  fam <- frmtmb_family(
    "rw_delta",
    dpars = c("alpha", "beta"),
    links = list(alpha = "logit", beta = "log"),
    primary_dpars = "alpha",
    type = "discrete",
    required_aterms = c("reward1", "reward2"),
    # The likelihood factorizes over trials but is not ROWWISE: a
    # trial's factor reads the whole history of its subject, which this
    # signature cannot see. Returning something anyway is the silent lie
    # the structured protocol exists to remove.
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A rw_delta() trial's probability depends on every earlier ",
           "trial of the same subject, so the family has no row-wise ",
           "log-density. Use logLik() for the total, or fitted() for ",
           "the per-trial choice probabilities", call. = FALSE)
    },
    valid_y = function(y, aterms) {
      if (!all(y %in% c(0, 1))) {
        stop("rw_delta(): the response is the arm chosen on each trial, ",
             "coded 1 for the first arm and 0 for the second",
             call. = FALSE)
      }
    },
    init_dpars = list(alpha = function(y, aterms) 0.3,
                      beta = function(y, aterms) 1),
    structure = rw_structure()
  )
  fam[["rw"]] <- list(subject_expr = subject_expr, trial_expr = trial_expr)
  fam
}
```

``` r

# The structure is the family's half of the protocol: where its data
# block comes from, what its likelihood is, what it can answer
# afterwards, and what it refuses. Every capability flag starts FALSE
# for a likelihood that does not factorize over rows, so the refusals
# below are explanations, not switches.
rw_structure <- function() {
  frmtmb_structure(
    frame_vars = function(fam) {
      list(fam[["rw"]][["subject_expr"]], fam[["rw"]][["trial_expr"]])
    },
    # An unanswered trial produces no prediction error, so dropping its
    # row IS the right recursion. This is where a delta rule differs
    # from a hidden Markov chain, which keeps the row because the state
    # still moves; see frmtmb.latent::hmm().
    keep_na = FALSE,
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("rw_delta() supports univariate models only: the recursion ",
             "is a likelihood over one response's trial sequences",
             call. = FALSE)
      }
      bad <- intersect(c("weights", "cens", "trunc_lb", "trunc_ub", "se"),
                       names(av))
      if (length(bad)) {
        stop("rw_delta() cannot be combined with ", bad[1L], "(): that ",
             "term reshapes a per-row likelihood contribution, and a ",
             "trial's contribution here is conditional on every earlier ",
             "trial of the same subject", call. = FALSE)
      }
    },
    frame_block = rw_block,
    loglik = rw_loglik,
    unit = "one subject's trial sequence",
    fitted_mean = rw_fitted,
    fitted_var = function(fit, block) {
      p <- rw_fitted(fit, block)
      p * (1 - p)
    },
    sim_ctx = rw_sim,
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a rw_delta() fit is not ",
        "available for newdata: a trial's choice probability is ",
        "conditional on the payoffs and choices of every earlier trial ",
        "of the same subject, and newdata carries no block to replay. ",
        "Use fitted() on the training data, or predict(dpar = 'alpha') ",
        "for a learning rate"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a rw_delta() fit: ",
        "the expected response is a choice probability that depends on ",
        "a whole trial history, which the synthetic grid this function ",
        "builds does not have. Plot the learning rate itself with ",
        "predict(dpar = 'alpha')"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a rw_delta() fit: ",
        "the tape holds the recursion over each whole sequence with no ",
        "registered observation vector. Use type = 'pearson', which ",
        "divides by the binomial variance of the choice probability"),
      deviance = paste0(
        "residuals(type = 'deviance') is not available for a rw_delta() ",
        "fit: the per-trial factors exist, but the structured ",
        "protocol's loglik slot returns one total, so the core never ",
        "sees them. Use type = 'pearson'")
    )
  )
}
```

The refusals are the interesting part. A structured family starts with
every capability turned off, and each one it leaves off should say why
in its own words.
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
is refused not because it is hard but because it is undefined: the
expected response here is a choice probability that depends on a whole
trial history, and the synthetic grid that function builds has no
history.

## Fitting it

``` r

# Scaffolding for the worked example, not part of the family. It builds
# a two-armed bandit with the payoff schedule of BOTH arms fixed in
# advance, one row per subject and trial, and a placeholder response the
# simulator overwrites. Arm 1 pays with probability p1 and arm 2 with
# probability p2, so the learnable contrast is p1 - p2.
rl_bandit_design <- function(n_subj, n_trial, p1 = 0.7, p2 = 0.3) {
  dd <- expand.grid(trial = seq_len(n_trial), id = seq_len(n_subj))
  dd$condition <- factor(rep(0:1, length.out = n_subj)[dd$id], 0:1,
                         c("ctl", "trt"))
  dd$id <- factor(dd$id)
  dd$pay1 <- stats::rbinom(nrow(dd), 1L, p1)
  dd$pay2 <- stats::rbinom(nrow(dd), 1L, p2)
  dd$choice <- 0
  # expand.grid() leaves an out.attrs attribute that str() prints at
  # length and nothing reads
  attr(dd, "out.attrs") <- NULL
  dd
}

# The learning rate carries the condition effect; the softmax
# temperature is a per-subject intercept correlated with it.
rl_bform <- bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
               beta ~ 1 + (1 | p | id))

rl_family <- function() rw_delta(subject = id, trial = trial)

# Parameter values in frm_simulate()'s natural spelling: the same names
# fixef() and VarCorr() report back.
rl_truth <- list(
  alpha_Intercept = stats::qlogis(0.35),
  alpha_conditiontrt = 0.8,
  beta_Intercept = log(3),
  sd_id__choice.alphaIntercept = 0.5,
  sd_id__choice.betaIntercept = 0.4,
  cor_id__choice.alphaIntercept__choice.betaIntercept = 0.3
)

# Draws come from the family's own simulator: frm_simulate() draws a new
# set of subject effects from `truth`, then walks each subject forward.
rl_simulate <- function(design, truth = rl_truth, nsim = 1, seed = NULL) {
  sims <- frm_simulate(rl_bform, data = design, family = rl_family(),
                       newparams = truth, nsim = nsim, seed = seed)
  lapply(seq_len(nsim), function(s) {
    d <- design
    d$choice <- sims[[s]]
    d
  })
}
```

``` r

set.seed(2024)
design <- rl_bandit_design(n_subj = 30, n_trial = 80)
dd <- rl_simulate(design, nsim = 1, seed = 2024)[[1]]
str(dd[1:3, ])
#> 'data.frame':    3 obs. of  6 variables:
#>  $ trial    : int  1 2 3
#>  $ id       : Factor w/ 30 levels "1","2","3","4",..: 1 1 1
#>  $ condition: Factor w/ 2 levels "ctl","trt": 1 1 1
#>  $ pay1     : int  0 1 1
#>  $ pay2     : int  0 0 1
#>  $ choice   : num  1 0 0
```

`rl_simulate()` goes through
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md),
which draws a fresh set of subject effects from `rl_truth` and then
calls the family’s own simulator. The choices in `dd` therefore come
from the same code path the fit is about to score.

Now the fit. This is an ordinary
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) call:

``` r

fit <- frm(bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
              beta ~ 1 + (1 | p | id)),
           family = rw_delta(subject = id, trial = trial),
           data = dd)
fixef(fit)
#> $alpha
#>  (Intercept) conditiontrt 
#>   -0.7919201    1.1220594 
#> 
#> $beta
#> (Intercept) 
#>    1.058637
```

The main right-hand side is the LEARNING RATE’s predictor, because the
family declared `primary_dpars = "alpha"`. The `| p |` in both terms
correlates each subject’s learning rate with its own temperature, which
is the usual hierarchical parameterization for these models.

``` r

VarCorr(fit)
#>   alpha: 1 | id + beta: 1 | id [ID] 
#>                      Name Std.Dev. choice.alpha:(Intercept)
#>  choice.alpha:(Intercept)  0.46280                         
#>   choice.beta:(Intercept)  0.37934                  -0.0422
```

Compare with the simulated values:

``` r

rbind(estimate = c(unlist(fixef(fit)),
                   sd_alpha = sqrt(VarCorr(fit)[[1]][1, 1]),
                   sd_beta = sqrt(VarCorr(fit)[[1]][2, 2])),
      truth = c(rl_truth$alpha_Intercept, rl_truth$alpha_conditiontrt,
                rl_truth$beta_Intercept,
                rl_truth$sd_id__choice.alphaIntercept,
                rl_truth$sd_id__choice.betaIntercept))
#>          alpha.(Intercept) alpha.conditiontrt beta.(Intercept)  sd_alpha
#> estimate        -0.7919201           1.122059         1.058637 0.4628031
#> truth           -0.6190392           0.800000         1.098612 0.5000000
#>            sd_beta
#> estimate 0.3793421
#> truth    0.4000000
```

This is one dataset, so read the differences as sampling error and not
as bias. The recovery study below is the systematic version.

Everything downstream of the fit works because the structure said it
does. The fitted values are per-trial choice probabilities:

``` r

p <- fitted(fit)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02198 0.59464 0.76165 0.72271 0.89176 0.99444
# the model's mean predicted choice rate against the observed one
c(predicted = mean(p), observed = mean(dd$choice))
#> predicted  observed 
#> 0.7227135 0.7104167
```

## Check 1: an independent reference for the recursion

The taped recursion is vectorized and masked, which is exactly the kind
of code that can be wrong in a way that still converges. So write the
model out again, the way it is stated, one subject and one trial at a
time, over quantities pulled back through the public accessors, and
compare.

``` r

u <- ranef(fit)[[1]]
sid <- match(as.character(dd$id), rownames(u))
xa <- model.matrix(~ condition, dd)
alpha <- plogis(as.vector(xa %*% unname(fixef(fit)$alpha)) + u[sid, 1])
beta <- exp(unname(fixef(fit)$beta)[[1]] + u[sid, 2])

ll <- 0
for (s in unique(sid)) {
  rows <- which(sid == s)
  rows <- rows[order(dd$trial[rows])]
  q <- c(0, 0)
  for (i in rows) {
    ll <- ll + dbinom(dd$choice[i], 1, plogis(beta[i] * (q[1] - q[2])),
                      log = TRUE)
    k <- if (dd$choice[i] == 1) 1 else 2
    pay <- if (k == 1) dd$pay1[i] else dd$pay2[i]
    q[k] <- q[k] + alpha[i] * (pay - q[k])
  }
}
c(reference = ll,
  from_fitted = sum(dbinom(dd$choice, 1, p, log = TRUE)),
  difference = ll - sum(dbinom(dd$choice, 1, p, log = TRUE)))
#>     reference   from_fitted    difference 
#> -1.149687e+03 -1.149687e+03  2.728484e-12
```

`tests/testthat/test-rl-example.R` runs this check on an unbalanced
design as well, where the padding is live, and against the full joint
density rather than the data part alone.

## Check 2: an identity against Stan

The reference above shares R with the thing it checks. A second
implementation in another language and another autodiff system does not.
`tests/testthat/helper-rl.R` holds a thirty-line Stan program for the
same model, written from the equations at the top of this page, and
compares Stan’s `log_prob` with frmtmb’s joint log density at frmtmb’s
own estimates. The Stan program declares exactly the parameters frmtmb
estimates, on the scales frmtmb estimates them on, so the map between
them is the identity and the comparison carries no Jacobian.

Measured on 30 subjects and 100 trials:

| quantity | at the estimates | at a displaced point |
|----|----|----|
| Stan `log_prob` | -1427.8453940 | -1460.3042080 |
| frmtmb joint log density | -1427.8453940 | -1460.3042080 |
| difference | -4.5e-13 | -2.5e-12 |
| largest gradient on the subject effects | 1.2e-14 | not asserted |
| largest gradient overall | 5.2 | not asserted |

The difference is zero and is expected to be exactly zero, not merely
small: core
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) carries no
default prior, the Stan program uses `target +=` with `_lpmf` and
`_lpdf` so it drops no normalizing constant, and every declared
parameter is unbounded so `adjust_transform = FALSE` has no Jacobian to
discard. The second column is a deliberately non-stationary point, so
the agreement cannot be an artifact of both sides sitting at an optimum.

The gradient row is the check that catches a wrong map. frmtmb puts the
subject effects at their conditional modes, so Stan’s gradient with
respect to them must vanish there. The gradient with respect to the
fixed effects is NOT zero at that point and is not asserted on: the
maximum-likelihood estimate makes the MARGINAL likelihood stationary,
not the joint one.

The test is gated the same way the brms comparison tier is, on
`FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN`, because it compiles a Stan
program.

## Check 3: parameter recovery

Simulate 100 datasets at known parameter values through the family’s own
simulator, fit each one, and ask whether the estimates land on the truth
and whether their intervals cover it. The design is 40 subjects with 100
trials each and a condition effect on the learning rate. The script is
described in `dev/rl-findings.md`.

| parameter            | truth  | bias   | Monte Carlo se | sd of estimates | coverage |
|----------------------|--------|--------|----------------|-----------------|----------|
| `alpha_(Intercept)`  | -0.619 | +0.009 | 0.018          | 0.178           | 0.94     |
| `alpha_conditiontrt` | 0.800  | -0.018 | 0.025          | 0.248           | 0.94     |
| `beta_(Intercept)`   | 1.099  | +0.002 | 0.007          | 0.070           | 0.96     |
| `log sd(alpha)`      | -0.693 | -0.211 | 0.057          | 0.573           | 0.93     |
| `log sd(beta)`       | -0.916 | -0.062 | 0.017          | 0.167           | 0.95     |

The three fixed effects are unbiased to within their Monte Carlo error
and their 95 percent Wald intervals cover at the nominal rate. The
variance components are the weak point, and the learning-rate one
especially: its estimate is biased downward and scattered widely, which
is the usual behavior of a variance component estimated from binary data
with a small number of levels. Coverage is over the replicates that
produced a finite interval, which for `log sd(alpha)` was 98 of the 100.

## The Laplace caveat, measured

frmtmb integrates the subject effects out with a Laplace approximation.
The approximation is a Gaussian fitted at the conditional mode, and it
is exact only when the conditional log-density really is quadratic. For
a binary response it is not, and the fewer trials a subject has, the
less quadratic it is.

The usual way to price that error in frmtmb is `frm(importance =)`,
which reweights draws from the Laplace Gaussian. It is REFUSED here:

``` r

frm(bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
       beta ~ 1 + (1 | p | id)),
    family = rw_delta(subject = id, trial = trial),
    data = dd, importance = 200)
#> Error:
#> ! `importance` cannot correct the 'rw_delta' family: it supplies its own log-likelihood, which does not factorize over rows, so a group's rows have no separable integrand to resample. This is the same restriction quadrature has. Use importance = 0
```

That refusal is honest but wider than it needs to be, and the closing
section returns to it.

What can still be measured is the consequence. Run the same recovery
study at 20 trials per subject instead of 100, and read the bias and the
coverage:

| parameter            | 100 trials, bias | coverage | 20 trials, bias | coverage |
|----------------------|------------------|----------|-----------------|----------|
| `alpha_(Intercept)`  | +0.009           | 0.94     | +0.064          | 0.95     |
| `alpha_conditiontrt` | -0.018           | 0.94     | -0.005          | 0.95     |
| `beta_(Intercept)`   | +0.002           | 0.96     | -0.000          | 0.92     |
| `log sd(alpha)`      | -0.211           | 0.93     | -0.628          | 0.73     |
| `log sd(beta)`       | -0.062           | 0.95     | -0.247          | 0.89     |

100 replicates each, 40 subjects, the same true parameter values.
Coverage is over the replicates that produced a finite interval: at 20
trials that was 84 of 100 for `log sd(alpha)` and 96 or more for
everything else.

Two readings. The fixed effects survive. At 20 trials per subject their
bias is still inside Monte Carlo error and their intervals still cover
close to the nominal rate, so a study that wants the condition effect on
the learning rate can work with short sessions.

The variance components do not survive. The learning rate’s log standard
deviation comes out 0.63 too low, which is a standard deviation about
half the true one, and its interval covers 73 times in 100. Sixteen of
the hundred short-session fits gave no usable interval for it at all.

Do not read all of that as Laplace error. A variance component estimated
by maximum likelihood from binary data with 40 levels is biased downward
whether or not the integral is approximated, and this table does not
separate the two causes. Separating them is what `frm(importance =)`
exists for, and this family cannot ask it. What the table does establish
is where the answers are safe: report fixed effects from short sessions,
and treat a subject-level standard deviation estimated from 20 binary
trials as a lower bound.

## Shape and cost

The family loops over TRIALS and vectorizes over subjects. It is worth
seeing what the other spelling costs. Here is the same likelihood
written the way the model is usually stated, one row at a time, with the
value store updated by sub-assignment:

``` r

for (t in seq_len(ncol(idx))) {
  for (s in seq_len(nrow(idx))) {          # the extra loop
    i <- idx[s, t]
    eta <- beta[i] * (q1[s] - q2[s])
    ll <- ll + y[i] * eta - RTMB::logspace_add(0 * eta, eta)
    if (y[i] == 1) q1[s] <- q1[s] + alpha[i] * (pay1[i] - q1[s])
    else q2[s] <- q2[s] + alpha[i] * (pay2[i] - q2[s])
  }
}
```

All three spellings give the same objective value to fourteen
significant digits. On 40 subjects and 100 trials, 4000 rows, with the
three builds interleaved and each gradient timed over a batch of 100
calls:

| `loglik` shape                        | tape build | one gradient |
|---------------------------------------|------------|--------------|
| vectorized over subjects              | 0.08 s     | 7.3 ms       |
| row loop, scalar value store          | 0.57 s     | 7.2 ms       |
| row loop, length-n log-density vector | 0.65 s     | 7.2 ms       |

The whole penalty is paid at TAPE CONSTRUCTION, and none of it per
evaluation. That is what it has to be: once the tape exists it is a node
list, and the R code that built it is gone. The gradient column is flat
to within measurement noise, and an optimization runs hundreds of
gradients against one build.

Do not read the build column as a fixed multiplier either. RTMB’s `[<-`
copies the vector it writes into, so the sub-assigning form gets
relatively worse as the response grows:

| rows  | vectorized | length-n sub-assignment | factor |
|-------|------------|-------------------------|--------|
| 1000  | 0.03 s     | 0.16 s                  | 5.3    |
| 4000  | 0.08 s     | 0.60 s                  | 7.5    |
| 10000 | 0.14 s     | 1.67 s                  | 11.9   |
| 20000 | 0.34 s     | 4.33 s                  | 12.7   |

Values agree to 1e-11 across the whole sweep. The rule to take away is
the shape, not a number: vectorize the tape-building loop, and expect
the cost of not doing so to grow with the data.

Scaling of the vectorized family, interleaved,
`frm(dry_run = "objective")` with the frame-assembly time subtracted,
median of 7 builds and the best of 3 gradient batches:

| subjects | trials | rows  | tape build | one gradient |
|----------|--------|-------|------------|--------------|
| 10       | 100    | 1000  | 0.04 s     | 10.0 ms      |
| 40       | 25     | 1000  | 0.03 s     | 8.8 ms       |
| 40       | 100    | 4000  | 0.08 s     | 39.6 ms      |
| 40       | 400    | 16000 | 0.28 s     | 153.4 ms     |
| 160      | 100    | 16000 | 0.22 s     | 160.6 ms     |
| 40       | 800    | 32000 | 0.61 s     | 296.2 ms     |

Both costs track the ROW count, which is the node count. Read the two
16000-row rows against each other for the part that does not: more
trials at a fixed row count costs more to tape, because the R loop runs
more iterations, and costs the same per gradient, because the tape is
the same size either way. The size of that build gap is not stable
enough to quote. Two interleaved runs of this table on one machine put
it at 1.27x and 1.70x.

## What a built-in family would add

The example above is a complete family. Three things would still change
if it moved into a package.

**A per-unit log-likelihood.** This family’s likelihood factorizes over
trials, but `frmtmb_structure(loglik =)` returns one AD scalar, so the
core never sees the individual factors. Two features follow from that
and are unavailable here: a pointwise log-likelihood matrix, which
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md) need, and
`frm(importance =)`, whose refusal appears above. The importance
correction needs one value per GROUP, not per row, and this family has
that. The refusal is written in terms of rows because that is the only
granularity the slot exposes. A future `pointwise_loglik` slot, or a
`loglik` allowed to return a vector over the structure’s own `unit`,
would serve both.

**Newdata.** `predict(type = "response")` on new data would need a block
built without a response, and the protocol defers that case on purpose.
For this family the meaning is clear enough (replay the new subject’s
payoff schedule against the fitted parameters), so a built-in version
would be worth the work.

**More of the literature.** A real package would carry the family of
models this one is the simplest member of: separate learning rates for
gains and losses, a decay term, more than two arms. All of them share
this block and this loop shape, and differ only in the body.

None of that is a limitation of the protocol as a place to put a
sequential likelihood. The core knows nothing about reinforcement
learning, and the model above fits, predicts, simulates and recovers its
parameters through the ordinary grammar.
