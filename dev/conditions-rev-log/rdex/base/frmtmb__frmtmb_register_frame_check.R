.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_register_frame_check
### Title: Check an assembled model frame from another package
### Aliases: frmtmb_register_frame_check

### ** Examples

# the shape of a check: refuse what only the assembled frame shows
check_not_constant <- function(spec, frame) {
  if (length(unique(frame$y)) == 1L) {
    stop("The response takes one value, so nothing can be estimated.",
         call. = FALSE)
  }
}

# a registration lasts for the session, so it belongs in the
# contributing package's .onLoad(), not in a script
## Not run: 
##D frmtmb_register_frame_check(check_not_constant)
## End(Not run)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
