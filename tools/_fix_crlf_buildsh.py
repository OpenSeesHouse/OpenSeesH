#!/usr/bin/env python3
from pathlib import Path
import sys
paths = sys.argv[1:] or ["/home/ubuntu/hpc1_build_rootcause.sh"]
for s in paths:
    p = Path(s)
    t = p.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    p.write_bytes(t)
    print("normalized", p)
