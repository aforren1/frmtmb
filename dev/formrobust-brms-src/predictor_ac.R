predictor_ac <- 
function (eta, prep, i, fprep = NULL) 
{
    if (!is.null(prep$ac[["err"]])) {
        eta <- eta + p(prep$ac$err, i, row = FALSE)
    }
    else if (has_ac_class(prep$ac$acframe, "arma")) {
        if (!is.null(i)) {
            stop2("Pointwise evaluation is not possible for ARMA models.")
        }
        eta <- .predictor_arma(eta, ar = prep$ac$ar, ma = prep$ac$ma, 
            Y = prep$ac$Y, J_lag = prep$ac$J_lag, fprep = fprep)
    }
    if (has_ac_class(prep$ac$acframe, "car")) {
        eta <- eta + .predictor_re(Z = p(prep$ac$Zcar, i), r = prep$ac$rcar)
    }
    eta
}
