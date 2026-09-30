scale_unit <- 
function (x, lb = min(x), ub = max(x)) 
{
    (x - lb)/(ub - lb)
}
