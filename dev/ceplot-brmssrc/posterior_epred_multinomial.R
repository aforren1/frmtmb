posterior_epred_multinomial <- 
function (prep) 
{
    get_counts <- function(i) {
        eta <- insert_refcat(slice_col(eta, i), refcat = prep$refcat)
        dcategorical(cats, eta = eta) * trials[i]
    }
    eta <- get_Mu(prep)
    cats <- seq_len(prep$data$ncat)
    trials <- prep$data$trials
    out <- abind(lapply(seq_cols(eta), get_counts), along = 3)
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- prep$cats
    out
}
