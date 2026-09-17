## Bernoulli message agreement between brms and the lane in the design run.
br <- readRDS("dev/famlink-rev-designs-brms.rds"); la <- readRDS("dev/famlink-rev-designs-lane.rds")
k <- names(la)
mb <- sapply(k, function(i) br[[i]]$n_bern_msg); ml <- sapply(k, function(i) la[[i]]$n_bern_msg)
sb <- sapply(k, function(i) br[[i]]$status); sl <- sapply(k, function(i) la[[i]]$status)
cat("brms messages:", sum(mb > 0), " of which lane fits:", sum(mb > 0 & sl == "ok"),
    " lane count equal there:", sum(mb > 0 & sl == "ok" & ml == mb), "\n")
cat("both fit, brms silent:", sum(mb == 0 & sl == "ok" & sb == "ok"),
    " lane messages there:", sum(ml[mb == 0 & sl == "ok" & sb == "ok"] > 0), "\n")
