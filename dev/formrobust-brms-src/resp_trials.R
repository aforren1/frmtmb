resp_trials <- 
function (x) 
{
    trials <- deparse0(substitute(x))
    class_resp_special("trials", call = match.call(), vars = nlist(trials))
}
