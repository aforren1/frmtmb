ma <- 
function (time = NA, gr = NA, q = 1, cov = FALSE) 
{
    label <- deparse0(match.call())
    time <- deparse0(substitute(time))
    gr <- deparse0(substitute(gr))
    .arma(time = time, gr = gr, p = 0, q = q, cov = cov, label = label)
}
