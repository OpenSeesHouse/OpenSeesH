#!/bin/bash
# Stage the tagged OpenSeesH .so into openseeshpy prebuilt/, build the local wheel, install it.
set -euo pipefail
export HOME=/home/ubuntu
export PATH="$HOME/.local/bin:$PATH"
LOG=/home/ubuntu/openseeshpy_wheel.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

echo "==== $(date) openseeshpy wheel start ===="
cd /home/ubuntu
rm -rf openseeshpy
mkdir -p openseeshpy
tar -xzf /home/ubuntu/openseeshpy_tree.tar.gz -C openseeshpy
cd /home/ubuntu/openseeshpy
# Tarball comes from a Windows clone, so shell scripts carry CRLF.
find . -name '*.sh' -exec sed -i 's/\r$//' {} +
echo "--- repo root:"
ls -la

SO=/home/ubuntu/dist/OpenSeesH_312.cpython-312-x86_64-linux-gnu.so
test -f "$SO"
mkdir -p prebuilt/linux_x86_64/cp312
cp -f "$SO" prebuilt/linux_x86_64/cp312/
echo "--- prebuilt:"
ls -lah prebuilt/linux_x86_64/cp312/

# Ubuntu 24.04 system Python is externally managed; keep tooling in user site.
# setuptools>=77 validates license expressions through packaging>=24.2
python3.12 -m pip install --user --break-system-packages --disable-pip-version-check -q -U \
  build wheel setuptools 'packaging>=24.2'
python3.12 -c 'import build, setuptools, wheel, packaging; print("BUILD_TOOLS_OK", setuptools.__version__, packaging.__version__)'

chmod +x scripts/build_wheels.sh
./scripts/build_wheels.sh

echo "--- dist:"
ls -la dist
WHL=$(ls dist/*cp312*linux_x86_64.whl | head -1)
python3.12 -m pip install --user --break-system-packages --force-reinstall "$WHL"
python3.12 - <<'PY'
import OpenSeesHpy
print("PKG_OK", getattr(OpenSeesHpy, "__version__", "n/a"), OpenSeesHpy.__file__)
PY
echo ALL_DONE
