.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-conditions
### Title: Classed conditions
### Aliases: frmtmb-conditions frm_stop frm_warning frm_message
###   frm_match_arg frm_family_package

### ** Examples

d <- data.frame(y = c(1, 2, 3))
e <- tryCatch(frm(y ~ 1, data = d, family = "not_a_family"),
              frmtmb_error = function(e) e)
class(e)
conditionMessage(e)

f <- function() frm_warning("a warning from f()")
w <- tryCatch(f(), frmtmb_warning = function(w) w)
class(w)
conditionCall(w)

g <- function(type = c("response", "link")) frm_match_arg(type)
g("resp")
e <- tryCatch(g("bogus"), frmtmb_error = function(e) e)
conditionMessage(e)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
