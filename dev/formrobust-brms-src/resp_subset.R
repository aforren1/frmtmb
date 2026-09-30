resp_subset <- 
function (x) 
{
    subset <- deparse0(substitute(x))
    class_resp_special("subset", call = match.call(), vars = nlist(subset))
}
