Adds drift-diffusion response time families to 'frmtmb', for two-choice
response time data. A noisy evidence accumulator runs between two
boundaries and the response time is the first time it touches one; the
drift rate, boundary separation, non-decision time and start point each
take their own regression formula, and Ratcliff's across-trial
variability in drift rate, start point and non-decision time is
available as three more. The parameterization follows 'brms'. The Wiener
family uses the Navarro and Fuss (2009) pair of series, written in plain
R so that it differentiates exactly on an 'RTMB' tape, with a smooth
weight between the two series because a tape cannot choose between them
on a parameter. The generalized family of Shinn, Lam and Murray (2020)
allows a state-dependent drift and moving boundaries, for which no
closed-form density exists; it solves the Fokker-Planck equation on a
grid that a change of variable holds fixed while the boundaries move,
which is what keeps the likelihood differentiable. The linear ballistic
accumulator of Brown and Heathcote (2008) is a race between any number
of accumulators, and reaches choices with more than two alternatives,
which a diffusion between two boundaries cannot.
