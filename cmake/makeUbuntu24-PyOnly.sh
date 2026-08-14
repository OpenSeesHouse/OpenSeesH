#!/bin/bash
# Build OpenSeesPy (opensees.so) on Ubuntu 24.04 for Python 3.12, 3.13, and 3.14.
# No Tcl. Uses conan-pyonly.txt (hdf5, zlib, eigen only).
#
# Prerequisites (example):
#   sudo apt-get install -y build-essential cmake ninja-build gfortran \
#     liblapack-dev libblas-dev python3-pip \
#     python3.12 python3.12-dev python3.13 python3.13-dev \
#     python3.14 python3.14-dev   # use deadsnakes if needed
#   pip3 install --user "conan<2"
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

CONAN_BIN="${CONAN_BIN:-$HOME/.local/bin/conan}"
JOBS="${JOBS:-$(nproc)}"
PY_VERSIONS="${PY_VERSIONS:-3.12 3.13 3.14}"

mkdir -p dist

for ver in $PY_VERSIONS; do
  build="build-py${ver}"
  py_exe="/usr/bin/python${ver}"
  py_inc="/usr/include/python${ver}"
  py_lib="/usr/lib/x86_64-linux-gnu/libpython${ver}.so"

  if [[ ! -x "$py_exe" ]]; then
    echo "SKIP python${ver}: $py_exe not found"
    continue
  fi
  if [[ ! -d "$py_inc" ]]; then
    echo "SKIP python${ver}: headers not found at $py_inc (install python${ver}-dev)"
    continue
  fi

  echo "==== Building OpenSeesPy for Python ${ver} ===="
  mkdir -p "$build"
  cd "$build"

  # Prefer py-only Conan recipe (no tcl): copy over conan.txt for this build tree
  if [[ -x "$CONAN_BIN" ]]; then
    if [[ -f "$ROOT/conan.txt" && ! -f "$ROOT/conan.txt.with_tcl" ]]; then
      cp -f "$ROOT/conan.txt" "$ROOT/conan.txt.with_tcl"
    fi
    cp -f "$ROOT/conan-pyonly.txt" "$ROOT/conan.txt"
    "$CONAN_BIN" install "$ROOT" --build=missing -if . || true
  fi

  cmake_args=(
    -S "$ROOT"
    -B .
    -DCMAKE_BUILD_TYPE=Release
    -DOPS_FINAL_TARGET=OpenSeesPy
    -DOPS_WITH_TCL=OFF
    -DPython_EXECUTABLE="$py_exe"
    -DPython_INCLUDE_DIRS="$py_inc"
  )
  if [[ -f "$py_lib" ]]; then
    cmake_args+=(-DPython_LIBRARIES="$py_lib")
  fi
  if [[ -f conanbuildinfo.cmake ]]; then
    :
  elif [[ -f conan_toolchain.cmake ]]; then
    cmake_args+=(-DCMAKE_TOOLCHAIN_FILE="$PWD/conan_toolchain.cmake")
  fi

  # Map version macro like Win64 OpenSeesPy configs
  case "$ver" in
    3.12) cmake_args+=(-DCMAKE_CXX_FLAGS="-DPy312") ;;
    3.13) cmake_args+=(-DCMAKE_CXX_FLAGS="-DPy313") ;;
    3.14) cmake_args+=(-DCMAKE_CXX_FLAGS="-DPy314") ;;
  esac

  cmake "${cmake_args[@]}"
  cmake --build . --config Release --target OpenSeesPy --parallel "$JOBS"

  outdir="$ROOT/dist/cp${ver//./}"
  mkdir -p "$outdir"
  if [[ -f OpenSeesPy.so ]]; then
    cp -f OpenSeesPy.so "$outdir/opensees.so"
  elif [[ -f lib/OpenSeesPy.so ]]; then
    cp -f lib/OpenSeesPy.so "$outdir/opensees.so"
  else
    so="$(find . -name 'OpenSeesPy.so' | head -n 1)"
    if [[ -z "$so" ]]; then
      echo "ERROR: OpenSeesPy.so not found for python${ver}"
      exit 1
    fi
    cp -f "$so" "$outdir/opensees.so"
  fi
  echo "Wrote $outdir/opensees.so"
  cd "$ROOT"
done

echo "Done. Artifacts under dist/cp312, dist/cp313, dist/cp314"
