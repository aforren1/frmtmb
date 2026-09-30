resp_cat <- 
function (x) 
{
    thres <- deparse0(substitute(x))
    str_add(thres) <- " - 1"
    class_resp_special("thres", call = match.call(), vars = nlist(thres, 
        gr = "NA"))
}
