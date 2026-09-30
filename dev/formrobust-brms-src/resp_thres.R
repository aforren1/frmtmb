resp_thres <- 
function (x, gr = NA) 
{
    thres <- deparse0(substitute(x))
    gr <- deparse0(substitute(gr))
    class_resp_special("thres", call = match.call(), vars = nlist(thres, 
        gr))
}
