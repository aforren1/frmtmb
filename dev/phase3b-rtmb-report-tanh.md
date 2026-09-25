# DRAFT, NOT FILED: the taped second derivative of tanh is NaN for |x| of 711 and more

Status: draft for the user to review. Do not file without their say so.
Address: TMB (github.com/kaskr/adcomp). Found through RTMB 2.0.

## Summary

With TMB 1.9.25 (through RTMB 2.0, R 4.6.1, Windows), the second
derivative of `tanh(x)` on a tape is NaN for |x| of 711 and more. The
true value is 0 to double precision. At 710 it is correct.

The cause is in `TMB/include/TMBad/global.hpp`, `struct TanhOp`, whose
`reverse()` is (line 3534 in the installed 1.9.25):

```cpp
args.dx(0) += args.dy(0) * (Type(1.) / (cosh(args.x(0)) * cosh(args.x(0))));
```

`cosh(x)` overflows to Inf past about 710.5. On the first derivative
that gives 1 / Inf = 0, which is right to double. On the second, the
tape differentiates this expression again, and the derivative of
1 / cosh^2 at an infinite cosh is Inf / Inf = NaN.

## Reproduction

```r
library(RTMB)
t4 <- MakeTape(function(x) tanh(x), 0)
h4 <- t4$jacfun()
for (x in c(700, 710, 711, 1e4, -711)) {
  cat(x, t4$jacobian(x), h4$jacobian(x), "\n")
}
```

Output (`dev/phase3b-log/rtmb-repro.txt`, from
`dev/phase3b-rtmb-repro.R`):

| x | d, tape | d2, tape |
|---|---|---|
| 700 | 0 | 0 |
| 710 | 0 | 0 |
| 711 | 0 | NaN |
| 1e4 | 0 | NaN |
| -711 | 0 | NaN |

True values: d = sech^2 x and d2 = -2 tanh x sech^2 x, both 0 to
double here.

## Why it matters

A logistic weight written as `0.5 * (1 + tanh(z))` saturates for large
z, and a model whose weight argument can reach 711 has a NaN Hessian,
so a Laplace approximation fails. frmtmb.eam clamps the argument at 40
(`ddm_tanh_s()` in `R/wiener-rtcdf.R`), where tanh is 1 to double.

## Suggested fix

Form the derivative from the output, which is bounded:

```cpp
args.dx(0) += args.dy(0) * (Type(1.) - args.y(0) * args.y(0));
```
