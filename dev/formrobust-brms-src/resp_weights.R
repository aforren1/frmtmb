resp_weights <- 
function (x, scale = FALSE) 
{
    weights <- deparse0(substitute(x))
    scale <- as_one_logical(scale)
    class_resp_special("weights", call = match.call(), vars = nlist(weights), 
        flags = nlist(scale))
}
