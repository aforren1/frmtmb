Runs on frmtmb 0.61.0 and later. On frmtmb 0.72.0 the standard errors
come from the exact Hessian. At a coherence of 1 that gives a finite,
very large standard error, where frmtmb 0.71.0 gave `NaN`. This version
refuses the interval that frmtmb.coupling 0.6.0 then returns with an
upper end of exactly 1, so update the two packages together.
