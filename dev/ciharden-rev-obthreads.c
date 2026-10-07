/* Reviewer probe: ask the libopenblas.dll this R process loaded how
   many threads it runs and how it was configured. */
#include <windows.h>
#include <R.h>
#include <Rinternals.h>

typedef int (*nt_fn)(void);
typedef char *(*cfg_fn)(void);

SEXP rev_ob_info(void) {
  HMODULE h = GetModuleHandleA("libopenblas.dll");
  SEXP out = PROTECT(allocVector(STRSXP, 3));
  if (!h) {
    SET_STRING_ELT(out, 0, mkChar("not loaded"));
    UNPROTECT(1);
    return out;
  }
  char path[MAX_PATH];
  GetModuleFileNameA(h, path, MAX_PATH);
  SET_STRING_ELT(out, 0, mkChar(path));
  nt_fn nt = (nt_fn) GetProcAddress(h, "openblas_get_num_threads");
  char buf[32];
  snprintf(buf, sizeof buf, "%d", nt ? nt() : -1);
  SET_STRING_ELT(out, 1, mkChar(buf));
  cfg_fn cf = (cfg_fn) GetProcAddress(h, "openblas_get_config");
  SET_STRING_ELT(out, 2, mkChar(cf ? cf() : "no config"));
  UNPROTECT(1);
  return out;
}
