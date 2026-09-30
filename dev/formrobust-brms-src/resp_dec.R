resp_dec <- 
function (x) 
{
    dec <- deparse0(substitute(x))
    class_resp_special("dec", call = match.call(), vars = nlist(dec))
}
