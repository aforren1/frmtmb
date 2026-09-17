.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-robust-dpars
### Title: Exact distributional parameters for a custom density
### Aliases: frmtmb-robust-dpars dpar_log dpar_log1m dpar_log_complement
###   dpar_complement

### ** Examples

# a gate at plogis(40), which is exactly 1 in double precision
on_tape <- list(zi = stats::plogis(40), .eta_zi = 40)
off_tape <- list(zi = stats::plogis(40))

# the value a density needs, and what the subtraction gives instead
dpar_log1m(on_tape, "zi", "logit")
log(1 - off_tape$zi)

# off the tape the entry is absent and the plain form comes back
dpar_log1m(off_tape, "zi", "logit")

# both terms from one log odds, for a density that needs the pair
unlist(dpar_log_complement(on_tape, "zi", "logit"))

# the log of a positive dpar, through that dpar's own link
dpar_log(list(shape = exp(3), .eta_shape = 3), "shape", "log")
dpar_log(list(shape = log1p(exp(3)), .eta_shape = 3), "shape",
         "softplus")

# a density written over the pair rather than over 1 - mu
d <- list(mu = stats::plogis(40), .eta_mu = 40)
mp <- dpar_complement(d, "mu", "logit")
c(p = mp$p, q = mp$q)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
