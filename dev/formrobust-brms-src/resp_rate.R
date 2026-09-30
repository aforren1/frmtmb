resp_rate <- 
function (denom) 
{
    denom <- deparse0(substitute(denom))
    class_resp_special("rate", call = match.call(), vars = nlist(denom))
}
