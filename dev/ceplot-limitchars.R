# Lane ceplot: brms 2.23.0's limit_chars() on the labels the hypothesis
# plot test uses, to pin frmtmb's hyp_limit_chars() against it.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
x <- c("(x+fb+z+x:fb+x:z)-(0.123456789) = 0", "(x) > 0",
       "(abcdefghijklmnopqrst) = 0", "(abcdefghijklmnopq) < 0")
dput(brms:::limit_chars(x, chars = 20))
dput(brms:::limit_chars(x, chars = 5))
