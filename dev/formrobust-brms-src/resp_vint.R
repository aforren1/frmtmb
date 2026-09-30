resp_vint <- 
function (...) 
{
    vars <- as.list(substitute(list(...)))[-1]
    class_resp_special("vint", call = match.call(), vars = vars)
}
