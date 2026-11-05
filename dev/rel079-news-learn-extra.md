
## Bug fixes

* `rlddm()` is a density, and declares `type = "continuous"`. A correct fit with
  tight response times, logLik +666.77, was told "a probability mass is at most
  1". DHARMa's `integerResponse` reads the same type, and is right on it now.
