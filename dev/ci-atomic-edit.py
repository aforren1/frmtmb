"""Lane ci: make every test write to the shared Stan cache atomic.
Defines stan_cache_write() in one helper per package and replaces each
saveRDS(<model>, path) that writes a compiled program. Each old text
must occur exactly once; line endings are kept as found."""
import pathlib

FUN = '''
# Write a compiled program to the Stan cache under its final name only
# once it is whole: saveRDS() to a temporary file in the same directory,
# then rename. A run killed during a plain saveRDS() left a truncated
# RDS under the final name, and CI now saves the cache of a job that
# failed or timed out, so such a file would reach later runs. The
# rename replaces an existing file on R 4.6 for Windows, but it can
# still fail there while another process holds the target open; the
# target is then a program another run wrote, so the copy in hand is
# dropped. Only a rename that fails with no target is an error.
stan_cache_write <- function(mod, path) {
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path),
                  fileext = ".part")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(mod, tmp)
  moved <- suppressWarnings(file.rename(tmp, path))
  if (!moved && !file.exists(path)) {
    stop("could not move ", tmp, " to ", path, call. = FALSE)
  }
  invisible(path)
}
'''


def edit(path, pairs, append_fun=False):
    p = pathlib.Path(path)
    b = p.read_bytes()
    crlf = b"\r\n" in b
    t = b.decode("utf-8").replace("\r\n", "\n")
    for old, new in pairs:
        assert t.count(old) == 1, (path, old)
        t = t.replace(old, new)
    if append_fun:
        assert "stan_cache_write <- function" not in t, path
        t = t.rstrip("\n") + "\n" + FUN
    if crlf:
        t = t.replace("\n", "\r\n")
    p.write_bytes(t.encode("utf-8"))
    print(path, "crlf" if crlf else "lf")


S = "extensions/frmtmb.sample/tests/testthat/"
edit("tests/testthat/helper-brms.R",
     [("  saveRDS(mod, path)\n", "  stan_cache_write(mod, path)\n")],
     append_fun=True)
edit("extensions/frmtmb.learn/tests/testthat/helper-stan.R",
     [("  saveRDS(mod, path)\n", "  stan_cache_write(mod, path)\n")],
     append_fun=True)
edit("extensions/frmtmb.spatial/tests/testthat/helper-brms.R",
     [("  saveRDS(mod, path)\n", "  stan_cache_write(mod, path)\n")],
     append_fun=True)
edit(S + "helper-brms.R", [], append_fun=True)
edit(S + "test-mo-coding.R",
     [("    saveRDS(mod, path)\n", "    stan_cache_write(mod, path)\n"),
      ("    saveRDS(mm, path)\n", "    stan_cache_write(mm, path)\n")])
edit(S + "test-loglik-autocor.R",
     [("      saveRDS(mod, path)\n", "      stan_cache_write(mod, path)\n")])
edit(S + "test-mi-draws-brms.R",
     [("    saveRDS(mod, path)\n", "    stan_cache_write(mod, path)\n")])

# test-default-priors-fit.R has the same line twice, at two sites
p = pathlib.Path(S + "test-default-priors-fit.R")
b = p.read_bytes()
crlf = b"\r\n" in b
t = b.decode("utf-8").replace("\r\n", "\n")
assert t.count("    saveRDS(mod, path)\n") == 2
t = t.replace("    saveRDS(mod, path)\n", "    stan_cache_write(mod, path)\n")
if crlf:
    t = t.replace("\n", "\r\n")
p.write_bytes(t.encode("utf-8"))
print(p, "crlf" if crlf else "lf", "(2 sites)")
