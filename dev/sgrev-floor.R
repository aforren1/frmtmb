# The version claim: frmtmb.sample as changed cannot LOAD against a
# released core, because .onLoad calls an export that core has not.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
B <- "C:/Users/adf44/source/r/rellib-r3"
F <- "C:/Users/adf44/source/r/sgrev-lib"
cat("released core version: ",
    as.character(packageVersion("frmtmb", lib.loc = B)), "\n")
cat("released core exports frm_install_generics: ",
    "frm_install_generics" %in%
      parseNamespaceFile("frmtmb", B)$exports, "\n")
cat("lane core exports it: ",
    "frm_install_generics" %in%
      parseNamespaceFile("frmtmb", F)$exports, "\n")
suppressMessages(loadNamespace("frmtmb", lib.loc = B))
cat("core loaded from: ", B, "\n")
r <- tryCatch({ suppressMessages(loadNamespace("frmtmb.sample",
       lib.loc = F)); "LOADED" },
     error = function(e) paste("FAILED:", conditionMessage(e)))
cat("lane frmtmb.sample against released core: ", r, "\n")
cat("DESCRIPTION floor in the worktree: ",
    grep("frmtmb", readLines(
      "extensions/frmtmb.sample/DESCRIPTION"), value = TRUE)[2], "\n")
