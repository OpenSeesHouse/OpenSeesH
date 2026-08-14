#!/bin/bash
set -euo pipefail
LOG=/home/ubuntu/openseesh_build_rootcause.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "==== ROOT CAUSE REPORT $(date -u) ===="
echo "Focus audit Built-Archived OK locally."
echo "OTHER libs present under HOME/lib."
echo "Compiler: -Wall -std=c++17 -Wno-overloaded-virtual, no -Werror."
echo "Sync: DAMPING_LIBS; no MATERIAL Py/ typo; BBar in ELE+ZHY; RENDERER has Renderer.o; BBar OPS_Stream*."
echo "==== $(date -u) incremental SRC+pythonmodule ===="
cd /home/ubuntu/OpenSeesH
export HOME=/home/ubuntu
make -C SRC/damping -j"$(nproc)"
make -C SRC/element/UP-ucsd -j"$(nproc)"
make -C SRC/material/uniaxial/PY -j"$(nproc)"
make -C SRC/api -j"$(nproc)"
make -C SRC/renderer -j"$(nproc)"
make -C SRC -j"$(nproc)"
echo "==== SRC done $(date -u) ===="
ls -lah /home/ubuntu/lib/libOpenSees.a
make -C SRC/interpreter pythonmodule -j"$(nproc)"
echo "==== interpreter done $(date -u) ===="
ls -lah SRC/interpreter/opensees.so
mkdir -p /home/ubuntu/dist
OUT=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
cp -f SRC/interpreter/opensees.so "$OUT"
ls -lah "$OUT"
python3.12 - <<'PY'
import importlib.util
p="/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so"
spec=importlib.util.spec_from_file_location("OpenSeesH_312", p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("IMPORT_OK")
PY
echo "==== ALL DONE $(date -u) ===="
