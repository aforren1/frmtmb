resp_trunc <- 
function (lb = -Inf, ub = Inf) 
{
    lb <- deparse0(substitute(lb))
    ub <- deparse0(substitute(ub))
    class_resp_special("trunc", call = match.call(), vars = nlist(lb, 
        ub))
}
