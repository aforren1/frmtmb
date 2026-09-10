# The code an `ndt_group()` column is keyed on

The addition-term registry coerces every `ndt_group()` column to a
numeric code before any family sees it, because an addition term's value
is baked into the tape as data. This is that coercion, exported so that
a caller holding a raw grouping can build the `aterms` list
[`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md)
reads without going through
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html).

## Usage

``` r
ndt_bound_key(x)
```

## Arguments

- x:

  A grouping: a factor, character, logical or numeric vector.

## Value

A numeric vector of codes, one per element of `x`.

## Details

The code is a function of the group's LABEL and of nothing else, so a
factor, a character vector, a logical and integer codes that name the
same groups all give the same bound, and none of them depends on the
level order of the frame the column was evaluated in.

## It is a hash, and it is not injective

The code is two polynomial rolling hashes over the label's code points,
packed into one double below `2^49`. Two different labels CAN therefore
share a code. Where both are in one column this refuses by name rather
than merging them. Where they are not, that is, where a label the fit
never saw collides with one it did, the unseen row is paired with the
collided group's bound instead of being refused, at a rate of about
1.8e-15 per pairing.

That residual is structurally zero for the identifier shapes anyone
uses: labels of one to three printable ASCII characters cannot collide
at all, and all 10^8 identifiers of the form `S00000000` hash to 10^8
distinct codes where a birthday count predicts 8.88 collisions.
Constructing an actual collision needed labels differing in four
positions over a span of 6401 code points.

## See also

[`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md),
which reads the codes under the name `ndt_group`.

## Examples

``` r
ndt_bound_key(factor(c("s1", "s2", "s1")))
#> [1] 253570877197 253587654414 253570877197
```
