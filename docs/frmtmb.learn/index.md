The value-learning models of the reinforcement-learning and
computational-psychiatry literature, written as 'frmtmb' families so
that each of their parameters is an ordinary distributional parameter
with its own linear predictor. A learning rate then takes a condition
effect, a smooth term or a correlated per-subject random effect the way
a mean does, and it comes back with a standard error, which is what a
separate model per group does not give. Six families share one recursion
that walks trials once and updates every subject at each step: the
two-armed delta learner and its dual-rate and counterfactual variants, a
Kalman filter over the restless four-armed bandit of Daw and others
(2006), the two-stage model-based and model-free hybrid of Daw and
others (2011), and a prospect-theory learner for the Iowa gambling task.
Families are named as 'hBayesDM' names them where a name exists. Each
one is checked by an identity against an independent 'Stan' program of
the same model at the same estimates, and by parameter recovery at a
realistic scale.
