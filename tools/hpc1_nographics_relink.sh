#!/bin/bash
set -euo pipefail
export HOME=/home/ubuntu
cd /home/ubuntu/OpenSeesH
rm -f SRC/recorder/OpsRecorderCommands.o SRC/interpreter/opensees.so
make -C SRC/interpreter pythonmodule -j"$(nproc)"
OUT=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
mkdir -p /home/ubuntu/dist
cp -f SRC/interpreter/opensees.so "$OUT"
echo "Non-Python undefined symbols:"
nm -D --undefined-only -C "$OUT" | grep -v '^ *U Py\|^ *U _Py\|^ *U Py_\|^ *w \|GLIBCXX\|GCC_\|CXXABI' || true
python3.12 - <<'PY'
import importlib.util
p="/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so"
spec=importlib.util.spec_from_file_location("OpenSeesH_312", p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("IMPORT_OK")
PY
ls -lah "$OUT"
echo ALL_DONE
