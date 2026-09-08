# One Stan program per family, written from the published equations.
#
# They share a preamble, because every family in this package shares a
# data layout: a subject-by-trial index matrix, a mask, a design matrix
# for the primary distributional parameter, and one random intercept on
# that parameter whose standard deviation rides along as data.
#
# The recursion below is written the way the model is STATED, one
# subject and one trial at a time with the value store updated by
# sub-assignment. That is deliberate: it is the spelling the package's
# own engine deliberately does not use, so an error in the vectorized,
# masked, indicator-selected form the engine does use cannot cancel.

ln_stan_head <- function(nd, extra = "", choice_lb = 1L) {
  paste(
    "data {",
    "  int<lower=1> N; int<lower=1> S; int<lower=1> T; int<lower=1> K;",
    "  array[S, T] int<lower=1> idx;",
    "  matrix[S, T] mask;",
    "  array[N] int<lower=1> subj;",
    # rlddm() codes its response column as the BOUNDARY reached, 0 or 1,
    # because that is what dec() reads; the six softmax families code it
    # as the option taken, 1 to K.
    paste0("  array[N] int<lower=", choice_lb, "> choice;"),
    "  matrix[N, K] X;",
    "  real<lower=0> sd_u;",
    extra,
    "}",
    "parameters {",
    paste0("  vector[K] b; vector[", nd, "] bd; vector[S] u;"),
    "}",
    "model {",
    "  vector[N] eta = X * b;",
    "  for (i in 1:N) eta[i] += u[subj[i]];",
    sep = "\n")
}

ln_stan_tail <- function() {
  paste("  target += normal_lpdf(u | 0, sd_u);", "}", sep = "\n")
}

ln_stan_code_delta <- function() {
  paste(
    ln_stan_head(1, "  vector[N] pay1; vector[N] pay2;"),
    "  vector[N] alpha = inv_logit(eta);",
    "  real tau = exp(bd[1]);",
    "  vector[S] q1 = rep_vector(0, S);",
    "  vector[S] q2 = rep_vector(0, S);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        vector[2] v;",
    "        v[1] = tau * q1[s]; v[2] = tau * q2[s];",
    "        target += v[choice[i]] - log_sum_exp(v);",
    "        if (choice[i] == 1) q1[s] += alpha[i] * (pay1[i] - q1[s]);",
    "        else q2[s] += alpha[i] * (pay2[i] - q2[s]);",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

ln_stan_code_dual <- function(split = "pe") {
  sel <- if (split == "outcome") {
    "        real w = (choice[i] == 1 ? pay1[i] : pay2[i]) > 0 ? 1.0 : 0.0;"
  } else {
    "        real w = pe > 0 ? 1.0 : 0.0;"
  }
  paste(
    ln_stan_head(2, "  vector[N] pay1; vector[N] pay2;"),
    "  vector[N] arew = inv_logit(eta);",
    "  real apun = inv_logit(bd[1]);",
    "  real tau = exp(bd[2]);",
    "  vector[S] q1 = rep_vector(0, S);",
    "  vector[S] q2 = rep_vector(0, S);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        vector[2] v;",
    "        v[1] = tau * q1[s]; v[2] = tau * q2[s];",
    "        target += v[choice[i]] - log_sum_exp(v);",
    "        real pe = choice[i] == 1 ? pay1[i] - q1[s] : pay2[i] - q2[s];",
    sel,
    "        real a = w * arew[i] + (1 - w) * apun;",
    "        if (choice[i] == 1) q1[s] += a * pe; else q2[s] += a * pe;",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

ln_stan_code_fict <- function() {
  paste(
    ln_stan_head(2, "  vector[N] pay1; vector[N] pay2;"),
    "  vector[N] alpha = inv_logit(eta);",
    "  real bias = bd[1];",
    "  real tau = exp(bd[2]);",
    "  vector[S] e1 = rep_vector(0, S);",
    "  vector[S] e2 = rep_vector(0, S);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        vector[2] v;",
    "        v[1] = tau * e1[s] + bias; v[2] = tau * e2[s];",
    "        target += v[choice[i]] - log_sum_exp(v);",
    "        real oc = choice[i] == 1 ? pay1[i] : pay2[i];",
    "        real t1 = choice[i] == 1 ? oc : -oc;",
    "        e1[s] += alpha[i] * (t1 - e1[s]);",
    "        e2[s] += alpha[i] * (-t1 - e2[s]);",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

ln_stan_code_kalman <- function(bonus = FALSE) {
  paste(
    ln_stan_head(if (bonus) 6 else 5, paste(
      "  matrix[N, 4] pay;",
      "  real<lower=0> sigma_o;", sep = "\n")),
    "  vector[N] tau = exp(eta);",
    if (bonus) "  real phi = bd[6];" else "",
    "  real lam = inv_logit(bd[1]);",
    "  real ctr = bd[2];",
    "  real mu0 = bd[3];",
    "  real sig0 = exp(bd[4]);",
    "  real sigd = exp(bd[5]);",
    "  matrix[S, 4] mu = rep_matrix(mu0, S, 4);",
    "  matrix[S, 4] vv = rep_matrix(sig0 ^ 2, S, 4);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        vector[4] v;",
    if (bonus) {
      paste0("        for (k in 1:4) v[k] = tau[i] * (mu[s, k]",
             " + phi * sqrt(vv[s, k]));")
    } else {
      "        for (k in 1:4) v[k] = tau[i] * mu[s, k];"
    },
    "        target += v[choice[i]] - log_sum_exp(v);",
    "        int c = choice[i];",
    "        real g = vv[s, c] / (vv[s, c] + sigma_o ^ 2);",
    "        mu[s, c] += g * (pay[i, c] - mu[s, c]);",
    "        vv[s, c] *= (1 - g);",
    "        for (k in 1:4) {",
    "          mu[s, k] = lam * mu[s, k] + (1 - lam) * ctr;",
    "          vv[s, k] = lam ^ 2 * vv[s, k] + sigd ^ 2;",
    "        }",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

ln_stan_code_pvl <- function() {
  paste(
    ln_stan_head(3, "  matrix[N, 4] pay;"),
    "  vector[N] alpha = inv_logit(eta);",
    "  real shape = inv_logit(bd[1]);",
    "  real lam = exp(bd[2]);",
    "  real tau = exp(bd[3]);",
    "  matrix[S, 4] q = rep_matrix(0, S, 4);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        vector[4] v;",
    "        for (k in 1:4) v[k] = tau * q[s, k];",
    "        target += v[choice[i]] - log_sum_exp(v);",
    "        int c = choice[i];",
    "        real x = pay[i, c];",
    "        real mag = abs(x) == 0 ? 0 : pow(abs(x), shape);",
    "        real ut = x >= 0 ? mag : -lam * mag;",
    "        q[s, c] += alpha[i] * (ut - q[s, c]);",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

ln_stan_code_ts <- function() {
  paste(
    ln_stan_head(6, paste(
      "  matrix[N, 4] pay;",
      "  array[N] int<lower=1, upper=2> s2;",
      "  array[N] int<lower=1, upper=2> c2;",
      "  real<lower=0.5, upper=1> pc;", sep = "\n")),
    "  vector[N] w = inv_logit(eta);",
    "  real a1 = inv_logit(bd[1]);",
    "  real t1 = exp(bd[2]);",
    "  real a2 = inv_logit(bd[3]);",
    "  real t2 = exp(bd[4]);",
    "  real lam = inv_logit(bd[5]);",
    "  real pers = bd[6];",
    "  matrix[S, 2] qmf = rep_matrix(0, S, 2);",
    "  matrix[S, 4] q2 = rep_matrix(0, S, 4);",
    "  matrix[S, 2] rp = rep_matrix(0, S, 2);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        real b1 = fmax(q2[s, 1], q2[s, 2]);",
    "        real b2 = fmax(q2[s, 3], q2[s, 4]);",
    "        vector[2] mb;",
    "        mb[1] = pc * b1 + (1 - pc) * b2;",
    "        mb[2] = (1 - pc) * b1 + pc * b2;",
    "        vector[2] v1;",
    "        for (k in 1:2)",
    "          v1[k] = t1 * (w[i] * mb[k] + (1 - w[i]) * qmf[s, k]",
    "                        + pers * rp[s, k]);",
    "        target += v1[choice[i]] - log_sum_exp(v1);",
    "        int base = 2 * (s2[i] - 1);",
    "        vector[2] v2;",
    "        for (k in 1:2) v2[k] = t2 * q2[s, base + k];",
    "        target += v2[c2[i]] - log_sum_exp(v2);",
    "        int m = base + c2[i];",
    "        real d1 = q2[s, m] - qmf[s, choice[i]];",
    "        real d2 = pay[i, m] - q2[s, m];",
    "        qmf[s, choice[i]] += a1 * (d1 + lam * d2);",
    "        q2[s, m] += a2 * d2;",
    "        rp[s, 1] = choice[i] == 1 ? 1 : 0;",
    "        rp[s, 2] = choice[i] == 2 ? 1 : 0;",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

# ORL. The two learning rates SWAP between the played deck and the three
# that were not, on the sign of what the played deck returned; getting
# that swap the wrong way round is the error this program exists to
# catch, so it is written as the published two-branch statement rather
# than as the arithmetic selector the engine uses.
ln_stan_code_orl <- function() {
  paste(
    ln_stan_head(4, "  matrix[N, 4] pay;"),
    "  vector[N] Arew = inv_logit(eta);",
    "  real Apun = inv_logit(bd[1]);",
    "  real kk = exp(bd[2]);",
    "  real betaF = bd[3];",
    "  real betaP = bd[4];",
    "  matrix[S, 4] ev = rep_matrix(0, S, 4);",
    "  matrix[S, 4] ef = rep_matrix(0, S, 4);",
    "  matrix[S, 4] ps = rep_matrix(0, S, 4);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        int c = choice[i];",
    "        vector[4] v;",
    "        for (k in 1:4)",
    "          v[k] = ev[s, k] + betaF * ef[s, k] + betaP * ps[s, k];",
    "        target += v[c] - log_sum_exp(v);",
    "        real x = pay[i, c];",
    "        real sg = x > 0 ? 1 : (x < 0 ? -1 : 0);",
    "        real a_play = x >= 0 ? Arew[i] : Apun;",
    "        real a_fic = x >= 0 ? Apun : Arew[i];",
    "        for (k in 1:4) {",
    "          if (k == c) {",
    "            ev[s, k] += a_play * (x - ev[s, k]);",
    "            ef[s, k] += a_play * (sg - ef[s, k]);",
    "            ps[s, k] = 1;",
    "          } else {",
    "            ef[s, k] += a_fic * (-sg / 3 - ef[s, k]);",
    "            ps[s, k] = ps[s, k] / exp(kk * log(3));",
    "          }",
    "        }",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}

# RLDDM. Stan has the Wiener first-passage density built in, which makes
# this the one program here that shares no code path with the R side at
# all: frmtmb.eam evaluates Navarro and Fuss's two series and blends
# them, and Stan uses its own implementation. So this row checks the
# DENSITY as well as the recursion, which none of the others do.
#
# Stan's wiener_lpdf is the UPPER-boundary density. The lower boundary
# is the same function at (1 - bias, -drift), which is the reflection
# the R side spells inside ddm_lpdf_both().
ln_stan_code_rlddm <- function() {
  paste(
    ln_stan_head(4, paste(
      "  vector[N] pay1; vector[N] pay2;",
      "  vector<lower=0>[N] rt;",
      "  real<lower=0> ndt_ub;", sep = "\n"), choice_lb = 0L),
    "  vector[N] alpha = inv_logit(eta);",
    "  real drift = bd[1];",
    "  real bs = exp(bd[2]);",
    "  real ndt = ndt_ub / (1 + exp(-bd[3]));",
    "  real bias = inv_logit(bd[4]);",
    "  vector[S] q1 = rep_vector(0, S);",
    "  vector[S] q2 = rep_vector(0, S);",
    "  for (t in 1:T) {",
    "    for (s in 1:S) {",
    "      if (mask[s, t] == 1) {",
    "        int i = idx[s, t];",
    "        real v = drift * (q2[s] - q1[s]);",
    "        if (choice[i] == 1) {",
    "          target += wiener_lpdf(rt[i] | bs, ndt, bias, v);",
    "          q2[s] += alpha[i] * (pay2[i] - q2[s]);",
    "        } else {",
    "          target += wiener_lpdf(rt[i] | bs, ndt, 1 - bias, -v);",
    "          q1[s] += alpha[i] * (pay1[i] - q1[s]);",
    "        }",
    "      }",
    "    }",
    "  }",
    ln_stan_tail(), sep = "\n")
}
