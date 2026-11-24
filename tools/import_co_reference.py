"""Import the authors' four coordinate blocks from Shimozono 2025 SI S31–S37.
Usage: python tools/import_co_reference.py local_si.pdf runs/co-reference
Requires Poppler pdftotext. Literature file stays local.
"""

import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


def extract(source, dest):
    text = subprocess.check_output(
        ["pdftotext", "-raw", "-f", "31", "-l", "38", str(source), "-"]
    ).decode()
    if "N-Phenyl Morpholine" not in text:
        raise ValueError("Not the expected SI / coordinate pages")
    lines = text.splitlines()
    results = {}
    pattern = r"^(H|C|N|B|O|Co)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)$"
    for label, expected in [("m-TS1", 96), ("p-TS1", 96), ("m-TS3", 89), ("p-TS3", 89)]:
        idx = lines.index(label)
        rows = []
        for line in lines[idx + 1 :]:
            line = line.strip()
            if not line or re.fullmatch(r"S\d+", line):
                continue
            if not re.fullmatch(pattern, line):
                break
            rows.append(line)
        comp = Counter(r.split()[0] for r in rows)
        if len(rows) != expected or comp["Co"] != 1 or comp["B"] != 1:
            raise ValueError(f"Unexpected coordinates for {label}: {len(rows)}, {comp}")
        results[label] = rows
    dest.mkdir(parents=True, exist_ok=False)
    for label, rows in results.items():
        (dest / (label + ".xyz")).write_text(
            f"{len(rows)}\nShimozono SI {label}; 10.1021/jacs.5c05965.s001\n"
            + "\n".join(rows)
            + "\n"
        )
    return results


if __name__ == "__main__":
    extract(Path(sys.argv[1]), Path(sys.argv[2]))
