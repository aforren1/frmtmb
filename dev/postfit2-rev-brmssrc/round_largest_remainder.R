function (x) 
{
    x <- as.numeric(x)
    total <- round(sum(x))
    out <- floor(x)
    diff <- x - out
    J <- order(diff, decreasing = TRUE)
    I <- seq_len(total - floor(sum(out)))
    out[J[I]] <- out[J[I]] + 1
    out
}
