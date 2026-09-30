arma <- 
function (time = NA, gr = NA, p = 1, q = 1, cov = FALSE) 
{
    label <- deparse0(match.call())
    time <- deparse0(substitute(time))
    gr <- deparse0(substitute(gr))
    .arma(time = time, gr = gr, p = p, q = q, cov = cov, label = label)
}
