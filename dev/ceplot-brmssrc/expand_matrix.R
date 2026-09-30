expand_matrix <- 
function (A, x, max_level = max(x), weights = 1) 
{
    stopifnot(is.matrix(A))
    stopifnot(length(x) == nrow(A))
    stopifnot(all(is_wholenumber(x) & x > 0))
    stopifnot(length(weights) %in% c(1, nrow(A), prod(dim(A))))
    A <- A * as.vector(weights)
    K <- ncol(A)
    i <- rep(seq_along(x), each = K)
    make_j <- function(n, K, x) K * (x[n] - 1) + 1:K
    j <- ulapply(seq_along(x), make_j, K = K, x = x)
    Matrix::sparseMatrix(i = i, j = j, x = as.vector(t(A)), dims = c(nrow(A), 
        ncol(A) * max_level))
}
