.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
suppressPackageStartupMessages(library(frmtmb))

a_function_with_a_long_name <- function(argument_one, argument_two) warning("a message that is long enough to wrap the prefix line")
a_function_with_a_long_name(argument_one = 1, argument_two = 2)
