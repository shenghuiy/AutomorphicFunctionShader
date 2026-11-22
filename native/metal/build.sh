#!/bin/sh
# Builds the native Metal renderers with Apple clang (Xcode or the Command Line Tools).
#
#   ./build.sh                       poincare, icosahedral, pearcey, kantorovich, inner and invariant_j in native/build
#   ./build.sh poincare              native/build/poincare only
#   ./build.sh poincare other        native/build/other (e.g. while poincare is running)
#
# The shaders are compiled at startup by the Metal framework, so the offline
# `metal` compiler is not needed. The CPU-side math headers are shared with the
# DirectX 12 programs in ../dx12.
set -e
cd "$(dirname "$0")"
OUT=../build
mkdir -p "$OUT"

build() {
    # $1 = program folder (also the source name), $2 = output name in $OUT
    echo "building $OUT/$2"
    clang++ -std=c++17 -O2 -Wall -fobjc-arc \
        -Icommon -I../dx12/common -I"../dx12/$1" \
        "$1/$1.mm" -o "$OUT/$2" \
        -framework Cocoa -framework Metal -framework QuartzCore
}

TARGET=${1:-all}
case "$TARGET" in
all)
    for p in poincare icosahedral pearcey kantorovich inner invariant_j; do build "$p" "$p"; done ;;
poincare|icosahedral|pearcey|kantorovich|inner|invariant_j)
    build "$TARGET" "${2:-$TARGET}" ;;
*)
    echo "Unknown target \"$TARGET\". Use poincare, icosahedral, pearcey, kantorovich, inner, invariant_j, or no argument for all." >&2
    exit 1 ;;
esac
