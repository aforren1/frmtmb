predictor.bprepl <- 
function (prep, i = NULL, fprep = NULL, ...) 
{
    nobs <- ifelse(!is.null(i), length(i), prep$nobs)
    eta <- matrix(0, nrow = prep$ndraws, ncol = nobs) + predictor_fe(prep, 
        i) + predictor_re(prep, i) + predictor_sp(prep, i) + 
        predictor_sm(prep, i) + predictor_gp(prep, i) + predictor_offset(prep, 
        i, nobs)
    eta <- predictor_ac(eta, prep, i, fprep = fprep)
    eta <- predictor_cs(eta, prep, i)
    unname(eta)
}
