mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/genrev-ab",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
say <- function(...) cat(sprintf(...))
if (mode == "ads-then-own") { q(library(genrevads)); q(loadNamespace("genrevown")) }
if (mode == "re-then-own")  { q(library(genrevre));  q(loadNamespace("genrevown")) }
if (mode == "ads-alone")    { q(library(genrevads)) }
if (mode == "own-then-ads") { q(library(genrevown)); q(library(genrevads)) }
if (mode == "own-reload")   { q(library(genrevads)); q(loadNamespace("genrevown"))
  q(unloadNamespace("genrevown")); q(loadNamespace("genrevown")) }
t <- structure(list(), class = "thing")
m <- structure(list(), class = "mine")
say("MODE %s\n", mode)
say("  generic owner      : %s\n",
    environmentName(topenv(environment(get("zap")))))
say("  zap(thing)         : %s\n",
    tryCatch(zap(t), error = function(e) paste("ERR", conditionMessage(e))))
say("  zap(mine)          : %s\n",
    tryCatch(zap(m), error = function(e) paste("ERR", conditionMessage(e))))
say("  getS3method thing  : %s\n",
    tryCatch(if (is.null(utils::getS3method("zap","thing",optional=TRUE)))
      "<none>" else "found", error = function(e) "ERR"))
say("  getS3method mine   : %s\n",
    tryCatch(if (is.null(utils::getS3method("zap","mine",optional=TRUE)))
      "<none>" else "found", error = function(e) "ERR"))
say("  binding active in ns: %s ; in package env: %s\n",
    bindingIsActive("zap", asNamespace("genrevads")),
    if ("package:genrevads" %in% search())
      bindingIsActive("zap", as.environment("package:genrevads")) else
    if ("package:genrevre" %in% search())
      paste("re:", bindingIsActive("zap", as.environment("package:genrevre")))
    else "n/a")
cat("GENREVDONE\n")
