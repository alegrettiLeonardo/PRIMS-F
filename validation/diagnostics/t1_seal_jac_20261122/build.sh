#!/bin/bash
# Build an offline diagnostic against the CURRENT worktree sources (release objects).
# Usage: build.sh PROGRAM.f90 OUTPUT_BINARY
set -e
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
make -C "$ROOT" -s -j4 release >/dev/null
API=$ROOT/validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90
B=$(mktemp -d)
OBJS=$(ls "$ROOT"/tmp/release/*.o | grep -v '/main.o$')
gfortran -std=f2018 -O2 -fopenmp -ffree-line-length-none -I"$ROOT/tmp/release" -J"$B" "$API" "$1" $OBJS -llapack -lblas -o "$2"
rm -rf "$B"
echo "built $2 from $(git -C "$ROOT" rev-parse --short HEAD)$(git -C "$ROOT" diff --quiet -- scr || echo +dirty)"
