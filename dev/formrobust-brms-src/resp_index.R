resp_index <- 
function (x) 
{
    index <- deparse0(substitute(x))
    class_resp_special("index", call = match.call(), vars = nlist(index))
}
