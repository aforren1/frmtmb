get_cond__ <- 
function (x) 
{
    out <- x[["cond__"]]
    if (is.null(out)) {
        out <- rownames(x)
    }
    as.character(out)
}
