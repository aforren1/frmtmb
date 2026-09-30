function (x) 
{
    x <- as.formula(x)
    if (length(x) == 3L) 
        update(x, . ~ 1)
    else NULL
}
