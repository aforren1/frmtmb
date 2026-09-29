# Construct the cases the guards in dev/release/build-docs.R exist to
# catch, and record what each one printed. Every guard this project has
# built failed OPEN on its first try, so a guard that has not been seen
# to fire is not evidence. Each case that is expected to FAIL is paired
# with a control that must PASS, because a check that fires on correct
# input is worse than no check.
#
# Run from anywhere. Writes dev/docsci-guard-tests.out through the
# caller's redirection, and one log per case under dev/.
$ErrorActionPreference = "Continue"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-docsci"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"

$env:R_LIBS = "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\rtools45\usr\bin;" + $env:PATH

function Case($name, $expect, $argv) {
  Write-Output ""
  Write-Output ("########## CASE " + $name + " (expect " + $expect + ")")
  $a = New-Object System.Collections.ArrayList
  [void]$a.Add("$ROOT/dev/release/build-docs.R")
  foreach ($x in $argv) { [void]$a.Add($x) }
  & $R $a.ToArray()
  $code = $LASTEXITCODE
  $got = if ($code -eq 0) { "pass" } else { "fail" }
  $verdict = if ($got -eq $expect) { "AS EXPECTED" } else { "WRONG" }
  Write-Output ("########## CASE " + $name + " exit " + $code +
                " -> " + $got + ", " + $verdict)
}

# --- C: a correct PARTIAL build. The control for B and D. One small
#     site, the minimum lowered to match. Nothing may fire.
Case "C-one-package-min-1" "pass" @(
  "--dest=$ROOT/dev/docsci-site-one",
  "--log=$ROOT/dev/docsci-one.log",
  "--pkgs=frmtmb.coupling",
  "--min-sites=1"
)

# --- D: the same one site with the minimum left at eight. The count arm
#     alone must fire, and it must say 1 of 8.
Case "D-one-package-min-8" "fail" @(
  "--dest=$ROOT/dev/docsci-site-one",
  "--log=$ROOT/dev/docsci-one-min8.log",
  "--pkgs=frmtmb.coupling"
)

# --- B: SEVEN packages, the default minimum of eight, which is the case
#     the task named. The seven extensions and not the core, so this
#     costs the cheap sites and not the expensive one.
Case "B-seven-packages-min-8" "fail" @(
  "--dest=$ROOT/dev/docsci-site-seven",
  "--log=$ROOT/dev/docsci-seven.log",
  "--pkgs=frmtmb.eam,frmtmb.latent,frmtmb.ode,frmtmb.spline,frmtmb.coupling,frmtmb.sample,frmtmb.learn"
)

# --- E: the clean refusal. A destination that holds files and is not a
#     pkgdown site must not be emptied.
$scratch = "$ROOT/dev/docsci-notasite"
Remove-Item -Recurse -Force $scratch -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $scratch | Out-Null
[System.IO.File]::WriteAllText("$scratch/precious.txt", "do not delete",
  (New-Object System.Text.UTF8Encoding($false)))
Case "E-clean-refusal" "fail" @(
  "--dest=$scratch",
  "--log=$ROOT/dev/docsci-refusal.log",
  "--pkgs=frmtmb.coupling",
  "--min-sites=1"
)
Write-Output ("E: precious.txt survived: " +
              (Test-Path "$scratch/precious.txt"))

# --- F: the article-gate guard, with a gate that CANNOT be satisfied.
#     A vignette is planted, because the script reads the gate names out
#     of the article sources and a planted name is the only honest way
#     to make one absent without touching a shared library. It goes in
#     the CORE vignettes directory and the cases build only
#     frmtmb.coupling, so the probe is SCANNED but never rendered: the
#     first spelling put it in the extension being built, and pkgdown
#     refused the file, which failed case G for a reason that had
#     nothing to do with the guard. It is removed again below, and the
#     removal is asserted.
$probe = "$ROOT/vignettes/docsci-gate-probe.Rmd"
[System.IO.File]::WriteAllText($probe,
  "---`ntitle: probe`n---`n`n``````{r}`nrequireNamespace(`"nosuch.docsci.pkg`")`n``````n",
  (New-Object System.Text.UTF8Encoding($false)))
# The destination is fresh, so the fact that it stays EMPTY is the
# evidence that the guard runs before the build and before the clean.
$gatedest = "$ROOT/dev/docsci-site-gate"
Remove-Item -Recurse -Force $gatedest -ErrorAction SilentlyContinue
Case "F-gate-absent-required" "fail" @(
  "--dest=$gatedest",
  "--log=$ROOT/dev/docsci-gate.log",
  "--pkgs=frmtmb.coupling",
  "--min-sites=1",
  "--require-articles"
)
$n = 0
if (Test-Path $gatedest) {
  $n = (Get-ChildItem -Force $gatedest | Measure-Object).Count
}
Write-Output ("F: entries left in the destination: " + $n + " (want 0)")

# --- G: the same planted gate WITHOUT --require-articles. It must warn
#     and build, because a local build should not be blocked by a
#     Suggests nobody has. This is F's control.
Case "G-gate-absent-not-required" "pass" @(
  "--dest=$ROOT/dev/docsci-site-gate2",
  "--log=$ROOT/dev/docsci-gate2.log",
  "--pkgs=frmtmb.coupling",
  "--min-sites=1"
)

Remove-Item -Force $probe
Write-Output ("probe removed: " + (-not (Test-Path $probe)))

# --- H: a destination whose last segment is not the package name. The
#     refusal happens during discovery, so no site is built and the case
#     is a second or two.
$yml = "$ROOT/extensions/frmtmb.coupling/_pkgdown.yml"
$orig = [System.IO.File]::ReadAllText($yml)
$bent = $orig.Replace("destination: ../../docs/frmtmb.coupling",
                      "destination: ../../docs/frmtmb.coupled")
if ($bent -eq $orig) { Write-Output "H: SKIPPED, could not bend the yml" }
else {
  [System.IO.File]::WriteAllText($yml, $bent,
    (New-Object System.Text.UTF8Encoding($false)))
  Case "H-destination-name-mismatch" "fail" @(
    "--dest=$ROOT/dev/docsci-site-bent",
    "--log=$ROOT/dev/docsci-bent.log",
    "--min-sites=1"
  )
  [System.IO.File]::WriteAllText($yml, $orig,
    (New-Object System.Text.UTF8Encoding($false)))
}
$restored = ([System.IO.File]::ReadAllText($yml) -eq $orig)
Write-Output ("H: _pkgdown.yml restored: " + $restored)

# --- I: a navbar href with no site behind it. Same method: bend the core
#     yml, run a FULL discovery against the eight-site tree that the main
#     run already built, with --no-clean so nothing is rebuilt.
$cyml = "$ROOT/_pkgdown.yml"
$corig = [System.IO.File]::ReadAllText($cyml)
$cbent = $corig.Replace("frmtmb/frmtmb.coupling/", "frmtmb/frmtmb.coupled/")
if ($cbent -eq $corig) { Write-Output "I: SKIPPED, could not bend the core yml" }
else {
  [System.IO.File]::WriteAllText($cyml, $cbent,
    (New-Object System.Text.UTF8Encoding($false)))
  Case "I-navbar-href-dead" "fail" @(
    "--dest=$ROOT/dev/docsci-site",
    "--log=$ROOT/dev/docsci-navbar.log",
    "--no-clean",
    "--pkgs=frmtmb.coupling",
    "--min-sites=1",
    "--check-links"
  )
  [System.IO.File]::WriteAllText($cyml, $corig,
    (New-Object System.Text.UTF8Encoding($false)))
}
Write-Output ("I: core _pkgdown.yml restored: " +
              ([System.IO.File]::ReadAllText($cyml) -eq $corig))
