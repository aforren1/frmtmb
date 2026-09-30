| outcome | 0.65.0 (rellib-r3) | lane build |
|---|---|---|
| pass | 275 | 306 |
| defect | 58 | 21 |
| divergence | 35 | 41 |
| cannot transfer | 126 | 126 |
| **total** | **494** | **494** |

Bin 1 passes: 275 of 494 before, 306 of 494 after.

Rows whose outcome or class moved: 38.

| row | before | after |
|---|---|---|
| `brmsfit-methods:747` | defect (argument) | defect (refuses-accepted) |
| `brmsfit-methods:283` | defect (argument) | divergence (class) |
| `brmsfit-methods:317` | defect (argument) | divergence (no-draws) |
| `brmsfit-methods:814` | defect (argument) | divergence (policy) |
| `brmsfit-methods:927` | defect (accepts-refused) | divergence (policy) |
| `data-helpers:12` | defect (different-error) | divergence (policy) |
| `priors:74` | defect (refuses-accepted) | divergence (policy) |
| `brm:102` | defect (internal-error) | pass |
| `brm:116` | defect (refuses-accepted) | pass (own-words) |
| `brm:81` | defect (different-error) | pass (own-words) |
| `brmsfit-methods:284` | defect (argument) | pass |
| `brmsfit-methods:285` | defect (output) | pass |
| `brmsfit-methods:301` | defect (refuses-accepted) | pass |
| `brmsfit-methods:305` | defect (argument) | pass |
| `brmsfit-methods:352` | defect (shape) | pass |
| `brmsfit-methods:737` | defect (refuses-accepted) | pass |
| `brmsfit-methods:741` | defect (refuses-accepted) | pass |
| `brmsfit-methods:772` | defect (argument) | pass |
| `brmsfit-methods:817` | defect (argument) | pass |
| `brmsfit-methods:82` | defect (argument) | pass |
| `brmsfit-methods:828` | defect (argument) | pass |
| `brmsfit-methods:83` | defect (argument) | pass |
| `brmsfit-methods:832` | defect (argument) | pass |
| `brmsfit-methods:840` | defect (refuses-accepted) | pass |
| `brmsfit-methods:841` | defect (refuses-accepted) | pass |
| `brmsfit-methods:904` | defect (output) | pass |
| `brmsfit-methods:924` | defect (output) | pass |
| `brmsfit-methods:938` | defect (internal-error) | pass (own-words) |
| `brmsfit-methods:94` | defect (argument) | pass |
| `brmsfit-methods:96` | defect (argument) | pass |
| `data-helpers:7` | defect (output) | pass |
| `families:102` | defect (argument) | pass |
| `priors:101` | defect (refuses-accepted) | pass |
| `priors:55` | defect (naming) | pass |
| `priors:59` | defect (naming) | pass |
| `standata:142` | defect (refuses-accepted) | pass |
| `standata:83` | defect (refuses-accepted) | pass |
| `standata:85` | defect (refuses-accepted) | pass |

Defect rows left: 21.

| row | class |
|---|---|
| `brmsfit-methods:154` | argument |
| `brmsfit-methods:205` | refuses-accepted |
| `brmsfit-methods:213` | refuses-accepted |
| `brmsfit-methods:215` | refuses-accepted |
| `brmsfit-methods:217` | output |
| `brmsfit-methods:269` | spelling |
| `brmsfit-methods:275` | spelling |
| `brmsfit-methods:277` | spelling |
| `brmsfit-methods:314` | argument |
| `brmsfit-methods:391` | argument |
| `brmsfit-methods:396` | argument |
| `brmsfit-methods:595` | argument |
| `brmsfit-methods:635` | naming |
| `brmsfit-methods:685` | filed |
| `brmsfit-methods:747` | refuses-accepted |
| `brmsfit-methods:953` | refuses-accepted |
| `brmsfit-methods:955` | refuses-accepted |
| `brmsfit-methods:959` | refuses-accepted |
| `brmsfit-methods:995` | spelling |
| `standata:75` | refuses-accepted |
| `standata:928` | internal-error |
