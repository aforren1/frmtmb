resp_vreal <- 
function (...) 
{
    vars <- as.list(substitute(list(...)))[-1]
    class_resp_special("vreal", call = match.call(), vars = vars)
}
