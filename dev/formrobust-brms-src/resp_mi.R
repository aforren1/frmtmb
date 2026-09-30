resp_mi <- 
function (sdy = NA) 
{
    sdy <- deparse0(substitute(sdy))
    class_resp_special("mi", call = match.call(), vars = nlist(sdy))
}
