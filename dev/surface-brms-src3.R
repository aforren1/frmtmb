# Lane surface: brms 2.23.0's categorical conditional-effects refusal.
source("dev/surface-env.R")
surface_env("base")
b <- deparse(body(get("conditional_effects.brmsfit", asNamespace("brms"))))
cat(grep("categorical", b, value = TRUE), sep = "\n")
