# The dev/test-backlog.md items lane illse closed, for dev/rel071-resolve.R
# (sourced only when illse is merged): the start of each item's first
# line, and the note that goes under it. dev/illse-findings.md has the
# numbers.
list(
  c("- **Accurate SEs on ill-conditioned designs** is its own lane: see",
    "CLOSED at 0.71.0 (lane illse): the covariance is computed in the design's coordinates and held as a factor, and every route frmtmb controls reads the factor; degree 8 with and without (1 | g) within 1.5e-5 of a 256-bit reference on the correlation scale (dev/illse-findings.md, sections 3 and 10)."),
  c("- **Accurate SEs on ill-conditioned designs.** A degree-8 raw",
    "CLOSED at 0.71.0 (lane illse), as above. marginaleffects and insight::get_varcov() get NaN with the frmtmb_cov_products warning where a product of the whole matrix loses its digits; vcov() keeps its entries and warns."),
  c("- **The convergence-warning policy**: the SE check does not run on a",
    "CLOSED at 0.71.0 (lane illse, item c): a stop nlminb did not call converged is judged as a converged fit when a Newton step from it gains at most grad_tol (conv_headroom()); otherwise its standard errors are kept with the reason stated. The xbeta fit of lane ci691 is judged on code 0 and code 1 alike."),
  c("- **The SE check does not run on a fit with a convergence warning**",
    "CLOSED at 0.71.0 (lane illse, item c), as above."),
  c("- **A residue of the exactly flat simplex direction off the",
    "CLOSED at 0.71.0 (lane illse, item b): se_mo_null() adds the split of an empty interior mo() category to the exact flat directions; 0 of 160 fits silent with finite zeta SEs (8 on 0.70.0; dev/illse-segap.R)."),
  c("- **Raw polynomials of degree 8 and up** are singular in double",
    "CLOSED at 0.71.0 (lane illse): degree 8 gets standard errors through the factor (lm()'s QR keeps it to 3.6e-9); degree 9 is where the frame drops an aliased column, so a degree-9 fit is a degree-8 one."),
  c("- **Gaussian linear predictors could take their SEs from a QR of the",
    "CLOSED at 0.71.0 (lane illse): reparam_basis() takes every dense design's QR, not only gaussian ones."),
  c("- **The optimizer stops short along an ill-conditioned direction**",
    "CLOSED at 0.71.0 (lane illse): newton_polish() steps in the design's coordinates where the natural Hessian needs them, also from a coded stop; the degree-8 fits end within a Newton gain of 1e-7 to 8e-7.")
)
