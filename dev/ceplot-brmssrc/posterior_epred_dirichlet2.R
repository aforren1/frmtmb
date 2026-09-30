posterior_epred_dirichlet2 <- 
function (prep) 
{
    mu <- get_Mu(prep)
    sums_mu <- apply(mu, 1:2, sum)
    cats <- seq_len(prep$data$ncat)
    for (i in cats) {
        mu[, , i] <- mu[, , i]/sums_mu
    }
    dimnames(mu)[[3]] <- prep$cats
    mu
}
