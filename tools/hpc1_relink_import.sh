#!/bin/bash
set -euo pipefail
export HOME=/home/ubuntu
cd /home/ubuntu/OpenSeesH
python3 - <<'PY'
from pathlib import Path
p = Path("Makefile.def")
t = p.read_text()
needle = "-D_LINUX -D_UNIX -D_PYTHON3 -D_CSS $(PY_VER_FLAG)"
repl = "-D_LINUX -D_UNIX -D_PYTHON3 -D_CSS -D_NO_NEW_RESTREPO $(PY_VER_FLAG)"
if "_NO_NEW_RESTREPO" in t:
    print("flags already have _NO_NEW_RESTREPO")
elif needle in t:
    p.write_text(t.replace(needle, repl, 1))
    print("patched _NO_NEW_RESTREPO")
else:
    raise SystemExit("flag pattern missing")
PY
rm -f SRC/interpreter/OpenSeesUniaxialMaterialCommands.o \
      SRC/interpreter/opensees.so \
      /home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
make -C SRC/interpreter pythonmodule -j"$(nproc)"
OUT=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
mkdir -p /home/ubuntu/dist
cp -f SRC/interpreter/opensees.so "$OUT"
ls -lah "$OUT"
python3.12 - <<'PY'
import importlib.util
p = "/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so"
spec = importlib.util.spec_from_file_location("OpenSeesH_312", p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("IMPORT_OK")
PY
if nm -C SRC/interpreter/opensees.so | grep -E 'U OPS_DoddRestr'; then
  echo "STILL_UNDEF_DODD"
  exit 1
fi
echo "ALL_DONE"
