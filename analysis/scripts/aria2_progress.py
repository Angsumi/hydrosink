#!/usr/bin/env python3
"""Print completed / total bytes from an aria2 .aria2 control file (format v1)."""
import struct, sys
b = open(sys.argv[1], "rb").read(); o = 6; ih = struct.unpack(">I", b[o:o+4])[0]; o += 4 + ih
pl = struct.unpack(">I", b[o:o+4])[0]; o += 4; tl = struct.unpack(">Q", b[o:o+8])[0]; o += 16
bl = struct.unpack(">I", b[o:o+4])[0]; o += 4
done = min(tl, sum(bin(x).count("1") for x in b[o:o+bl]) * pl)
print(f"{done/1e9:.2f}/{tl/1e9:.2f} GB ({100*done/tl:.0f}%)")
