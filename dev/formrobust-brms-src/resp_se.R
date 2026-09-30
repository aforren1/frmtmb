resp_se <- 
function (x, sigma = FALSE) 
{
    se <- deparse0(substitute(x))
    sigma <- as_one_logical(sigma)
    class_resp_special("se", call = match.call(), vars = nlist(se), 
        flags = nlist(sigma))
}
