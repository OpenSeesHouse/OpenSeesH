#!/bin/bash
# Smoke test the installed openseeshpy wheel (loads the tagged OpenSeesH .so).
set -euo pipefail
export HOME=/home/ubuntu
cd /home/ubuntu
python3.12 - <<'PY'
from OpenSeesHpy import OpenSeesHpy as oph

oph.wipe()
oph.model('basic', '-ndm', 2, '-ndf', 2)
oph.node(1, 0.0, 0.0)
oph.node(2, 100.0, 0.0)
oph.fix(1, 1, 1)
oph.fix(2, 0, 1)
oph.uniaxialMaterial('Elastic', 1, 2000.0)
oph.element('truss', 1, 1, 2, 10.0, 1)
oph.timeSeries('Linear', 1)
oph.pattern('Plain', 1, 1)
oph.load(2, 100.0, 0.0)
oph.system('BandGeneral')
oph.numberer('Plain')
oph.constraints('Plain')
oph.integrator('LoadControl', 1.0)
oph.algorithm('Linear')
oph.analysis('Static')
assert oph.analyze(1) == 0
d = oph.nodeDisp(2, 1)
print("nodeDisp(2,1) =", d)
expected = 100.0 * 100.0 / (2000.0 * 10.0)
assert abs(d - expected) < 1e-9, (d, expected)
print("SMOKE_OK")
PY
