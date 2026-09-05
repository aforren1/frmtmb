Two things a fitted curve needs after it is fitted. First, curve
inference: a smooth evaluated on a grid with pointwise and SIMULTANEOUS
confidence bands, its first and second derivatives with delta-method
standard errors, and the features movement papers report, the time of a
peak and the time a curve crosses a level, each with a standard error
from the implicit-function delta method. The simultaneous band is the
max-deviation simulation of Ruppert, Wand and Carroll (2003), the same
construction 'gratia' uses on 'mgcv' fits, and it is checked against
'gratia' inside its Monte Carlo error. Second, a Royston and Parmar
(2002) flexible parametric survival family, which writes the log
cumulative hazard as a natural cubic spline in log time. It is
parameterized exactly as 'flexsurv::flexsurvspline' parameterizes it, so
the two log likelihoods are the same number rather than two numbers that
ought to agree, and covariates reach the spline coefficients as
proportional hazards on the first and as time-varying effects on the
rest.
