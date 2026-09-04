# Compiling hmm.R and lca.R out of tree

Written before the files moved, which is the point of it. Both were
written inside the core package, where every internal helper is one
name away. An extension reaches only the exported surface, so each
internal name they touch had to be resolved first: already exported,
newly exported, or restructured away.

## Method

`codetools::findGlobals()` over every function whose `srcref` names
`R/hmm.R` or `R/lca.R` in a loaded core namespace, minus the two files'
own definitions, intersected with the core namespace and split on
`getNamespaceExports()`. The two files use no `:::`, no
`getFromNamespace()`, and no `do.call()` on a literal internal name, so
the static scan sees everything: verified by grep alongside the scan.

## Already exported: the step 9 accessors, and they were enough

`as_frmtmb_family()`, `compat_rule_builder()`, `dpar_linpred()`,
`eval_dpars()`, `frmtmb_family()`, `frmtmb_structure()`,
`latent_probs()`, `response_mean()`, `single_response()`.

Step 9 named this list from the protocol doc and exported it. The scan
confirms it, with one correction to the protocol's guesses:
`fit_extras()` and `linpred_key()` are NOT used by either file, and
`assemble_frame()` and `has_mixture()` are gone. `lca_profiles()` reads
the fit's own frame, as the protocol asked it to.

## Newly exported: four, each documented on `frmtmb-extension-api`

| symbol | who needs it | why not restructure |
|---|---|---|
| `frame_block_of(frame, resp)` | `hmm_parts()` | The read half of `frmtmb_structure(frame_block =)`. A family that writes a block at frame time must read it back at post-fit time, and `frame$blocks[[resp]]` is a layout an extension must not hardcode. Unavoidable for ANY structured family, not just these two. |
| `structure_supports_all(...)` | `lca_structure()` | The starting point for the second kind of structure the protocol describes: a family whose likelihood IS rowwise and which carries a structure only to hold two or three refusals. The all-`FALSE` default is exactly wrong for it. Spelling the TRUEs by hand would desync the moment core adds a capability flag. |
| `mixture_posterior(fit)` | `lca_structure()` | The protocol's stated reason `fam$mix` stayed on `lca()`: one posterior implementation serves `mixture()`, `mixture_mvn()` and `lca()`. Reimplementing it in the extension would fork the component interface's only consumer. |
| `mixture_multimodal_refusals(what)` | `lca_structure()` | Two user-facing sentences about a mathematical fact that is equally true of `mixture()` in core and `lca()` here. Copied text across a package boundary drifts silently; a core rewording must reach both. |

Nothing else was exported. In particular no `fam_structure()`: an
extension building a family already holds it, and reading another
family's is core's business.

## Restructured: three private clones, no new core surface

`%||%`, `check_flag()` and its `arg_desc()` are argument-plumbing
utilities, not extension API. Exporting a general TRUE/FALSE validator
would put it on the public surface forever to serve one call
(`lca(na.rm =)`). They are cloned verbatim into `R/utils.R` here
instead, so the message a user sees for `lca(na.rm = "yes")` is
unchanged to the character.

This is the message-uniqueness property working as designed: it holds
per package, not across the pair, so a template repeated here resolves
unambiguously to this source tree.

## Registrations at load

`hmm_compat_rules()` and `lca_compat_rules()` go to
`frmtmb::frmtmb_register_compat()` from `.onLoad()`.

No frame check is registered. `check_lca_structure()` runs through the
family's own `frmtmb_structure(check_frame =)` slot, which core calls
for whatever family the model names, so it needs no registry entry:
verified by reading `lca_structure()` and core's `check_frame` call
site. `frmtmb_register_frame_check()` is for a check that must run on
models the extension's family does NOT appear in, which is the ODE
package's case, not this one.
