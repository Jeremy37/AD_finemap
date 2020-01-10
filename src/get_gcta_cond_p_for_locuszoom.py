#!/usr/bin/env python
import sys
import re

for line in sys.stdin:
    vals = line.strip().split()
    
    marker = vals[1]
    if re.search("^rs", marker):
        marker = marker.split("_")[0]
    else:
        marker = vals[0] + ":" + vals[2]
    print marker + "\t" + vals[12]

