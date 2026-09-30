get_ac_vars <- 
function (x, var, ...) 
{
    var <- match.arg(var, c("time", "gr", "M"))
    acframe <- subset2(frame_ac(x), ...)
    out <- unique(acframe[[var]])
    setdiff(na.omit(out), "NA")
}
