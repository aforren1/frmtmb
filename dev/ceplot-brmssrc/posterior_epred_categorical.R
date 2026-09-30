posterior_epred_categorical <- 
function (prep) 
{
    get_probs <- function(i) {
        eta <- insert_refcat(slice_col(eta, i), refcat = prep$refcat)
        dcategorical(cats, eta = eta)
    }
    eta <- get_Mu(prep)
    cats <- seq_len(prep$data$ncat)
    out <- abind(lapply(seq_cols(eta), get_probs), along = 3)
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- prep$cats
    out
}
