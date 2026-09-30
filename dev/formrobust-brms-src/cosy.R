cosy <- 
function (time = NA, gr = NA) 
{
    label <- deparse0(match.call())
    time <- deparse0(substitute(time))
    time <- as_one_variable(time)
    gr <- deparse0(substitute(gr))
    stopif_illegal_group(gr)
    out <- nlist(time, gr, label)
    class(out) <- c("cosy_term", "ac_term")
    out
}
