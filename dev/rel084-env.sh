# Rscript under Git Bash with Rtools first on PATH segfaults
# (dev/round-handoff.md), so Rtools goes last. Sourced by the round's scripts.
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export GIT_OPTIONAL_LOCKS=0
