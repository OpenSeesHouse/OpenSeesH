#!/bin/bash
set -euo pipefail
export HOME=/home/ubuntu
cd /home/ubuntu/OpenSeesH
make -C SRC/element/surfaceLoad -j"$(nproc)"
# re-archive FE library with new objs (SRC all deletes then rearchives)
make -C SRC -j"$(nproc)"
rm -f SRC/interpreter/opensees.so
make -C SRC/interpreter pythonmodule -j"$(nproc)"
OUT=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
mkdir -p /home/ubuntu/dist
cp -f SRC/interpreter/opensees.so "$OUT"
ls -lah "$OUT"
echo "Undefined OPS_ symbols:"
nm -D --undefined-only -C "$OUT" | grep OPS_ || echo "(none)"
python3.12 - <<'PY'
import importlib.util
p="/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so"
spec=importlib.util.spec_from_file_location("OpenSeesH_312", p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("IMPORT_OK")
PY
echo ALL_DONE
