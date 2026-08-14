#!/bin/bash
# Root cause: Win64 links sdmuc.lib (SDM-UC.f); Linux never built/archived SDM-UC.o → undefined sdmuc_
set -euo pipefail
export HOME=/home/ubuntu
cd /home/ubuntu/OpenSeesH
LOG=/home/ubuntu/openseesh_sdmuc.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

echo "==== $(date) sdmuc root-cause relink ===="
test -f SRC/material/nD/stressDensityModel/SDM-UC.f
grep -q 'SDM-UC.o' SRC/Makefile
grep -q 'SDM-UC.o' SRC/material/nD/stressDensityModel/Makefile

rm -f SRC/material/nD/stressDensityModel/SDM-UC.o SRC/interpreter/opensees.so
make -C SRC/material/nD/stressDensityModel SDM-UC.o
ls -lah SRC/material/nD/stressDensityModel/SDM-UC.o
nm -g SRC/material/nD/stressDensityModel/SDM-UC.o | grep -i sdmuc || {
  echo "ERROR: sdmuc_ not in SDM-UC.o"
  exit 1
}

# Archive into FE library (same path as MATERIAL_LIBS AR)
ar r "$(grep '^FE_LIBRARY' Makefile.def | head -1 | sed 's/.*= *//')" \
  SRC/material/nD/stressDensityModel/SDM-UC.o 2>/dev/null \
  || ar r /home/ubuntu/lib/libOpenSees.a SRC/material/nD/stressDensityModel/SDM-UC.o
nm /home/ubuntu/lib/libOpenSees.a 2>/dev/null | grep -i 'sdmuc' | head -5 || true

make -C SRC/interpreter pythonmodule -j"$(nproc)"
OUT=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
mkdir -p /home/ubuntu/dist
cp -f SRC/interpreter/opensees.so "$OUT"
echo "Non-Python undefined (sample):"
nm -D --undefined-only -C "$OUT" | grep -vE '^ *U (Py|_Py|Py_)|^ *w |GLIBCXX|GCC_|CXXABI' | head -40 || true
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
