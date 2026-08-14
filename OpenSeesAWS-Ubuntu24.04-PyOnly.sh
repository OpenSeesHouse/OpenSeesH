#!/bin/bash
# Ubuntu 24.04 apt + CMake OpenSeesPy-only matrix (Python 3.12–3.14, no Tcl).
# Forked from OpenSeesAWS-Ubuntu22.04.sh style; does not build OpenSees/MP/SP.
set -euo pipefail

sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build gfortran \
  liblapack-dev libblas-dev python3-pip \
  python3.12 python3.12-dev \
  python3.13 python3.13-dev || true

# Python 3.14 may need deadsnakes on Ubuntu 24.04:
#   sudo add-apt-repository ppa:deadsnakes/ppa && sudo apt-get update
#   sudo apt-get install -y python3.14 python3.14-dev

pip3 install --user "conan<2" || true

ROOT="$(cd "$(dirname "$0")" && pwd)"
bash "$ROOT/cmake/makeUbuntu24-PyOnly.sh"
