# Package index

## Families

One recursion engine, eight models. Each parameter is an ordinary
distributional parameter with its own linear predictor, so a learning
rate takes a condition effect, a smooth term or a correlated per-subject
random effect the way a mean does. Seven turn a value store into a
softmax over options; rlddm() turns it into the drift rate of a
diffusion, so choices and response times are one likelihood.

- [`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)
  : Rescorla-Wagner delta learning on a two-armed bandit
- [`bandit2arm_dual()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_dual.md)
  : Two-armed delta learning with separate rates for gains and losses
- [`prl_fictitious()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md)
  : Counterfactual (fictitious) updating on two options
- [`bandit4arm2_kalman_filter()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit4arm2_kalman_filter.md)
  : A Kalman filter over the restless four-armed bandit
- [`ts_par7()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md)
  : The two-step task: a model-based and model-free hybrid
- [`igt_pvl_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_pvl_delta.md)
  : Prospect-valence learning with a delta rule, for the Iowa gambling
  task
- [`igt_orl()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_orl.md)
  : Outcome-representation learning for the Iowa gambling task
- [`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md)
  : Reinforcement learning with a drift-diffusion choice rule

## After the fit

These families declare no mean on the response scale, so fitted() and
residuals() refuse: for seven of them because the response is a nominal
option code, and for rlddm() because the mean of a first-passage time
belongs to frmtmb.eam. The trajectory is read here instead.

- [`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
  : Per-trial value estimates, prediction errors and choice
  probabilities
- [`frm_learn_families()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_learn_families.md)
  : The families this package supplies, and what each one is

## Simulation

Task designs with the payoff schedule of every option fixed in advance,
and draws from the generative process.

- [`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md)
  : Build a trial-level design for a learning task
- [`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
  : Draw datasets from a learning family's generative process

## Package

- [`frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
  [`frmtmb.learn-package`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
  : frmtmb.learn: Reinforcement-Learning Families for 'frmtmb' Models
