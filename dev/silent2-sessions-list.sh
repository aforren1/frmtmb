#!/usr/bin/env bash
# Lane silent2: run dev/silent2-session1.sh for each "<package> <order>"
# pair given, one after another.
#   bash dev/silent2-sessions-list.sh pkg1 order1 pkg2 order2 ...
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
while [ $# -ge 2 ]; do
  bash "$ROOT/dev/silent2-session1.sh" "$1" "$2"
  shift 2
done
