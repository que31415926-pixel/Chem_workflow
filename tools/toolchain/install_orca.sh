#!/usr/bin/env bash
# Verify (and on Linux x86-64, install) the locked ORCA 6.1.0 installer.
#
#   bash tools/toolchain/install_orca.sh /path/to/orca_6_1_0_linux_x86-64_shared_openmpi418.run [PREFIX]
#
# ORCA cannot be downloaded by script (registration + licence), so download the file named in
# toolchain.lock.json by hand from https://orcaforum.kofo.mpg.de first. This script then:
#   1. checks the file name, byte size and SHA-256 against toolchain.lock.json (refuses on mismatch);
#   2. on Linux x86-64: runs the installer into PREFIX (default $HOME/orca_6_1_0) and checks the
#      installed orca binary SHA-256 against the lock;
#   3. on other platforms: only verifies, and tells you to install that platform's 6.1.0 build;
#   4. writes PREFIX/chemflow_orca_receipt.json.
# Afterwards run: python3 tools/toolchain/check_toolchain.py --orca PREFIX/orca
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOCK="$REPO/toolchain.lock.json"
INSTALLER="${1:?usage: install_orca.sh /path/to/installer.run [PREFIX]}"
PREFIX="${2:-$HOME/orca_6_1_0}"

field() { python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); p=d["orca"]["platforms"]["linux-x86_64"]; print(p[sys.argv[2]])' "$LOCK" "$1"; }
sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
fsize() { if stat -c %s "$1" >/dev/null 2>&1; then stat -c %s "$1"; else stat -f %z "$1"; fi; }

WANT_NAME="$(field installer)"; WANT_BYTES="$(field installer_bytes)"
WANT_SHA="$(field installer_sha256)"; WANT_BIN="$(field orca_binary_sha256)"

[ -f "$INSTALLER" ] || { echo "FAIL installer not found: $INSTALLER"; exit 1; }
NAME="$(basename "$INSTALLER")"
[ "$NAME" = "$WANT_NAME" ] || { echo "FAIL file name $NAME != locked $WANT_NAME"; exit 1; }
BYTES="$(fsize "$INSTALLER")"
[ "$BYTES" = "$WANT_BYTES" ] || { echo "FAIL size $BYTES != locked $WANT_BYTES (incomplete download?)"; exit 1; }
echo "...  hashing $NAME"
SHA="$(sha256 "$INSTALLER")"
[ "$SHA" = "$WANT_SHA" ] || { echo "FAIL sha256 $SHA != locked $WANT_SHA"; exit 1; }
echo "PASS installer matches toolchain.lock.json ($NAME, $BYTES bytes)"

OS="$(uname -s)"; ARCH="$(uname -m)"
if [ "$OS" != "Linux" ] || [ "$ARCH" != "x86_64" ]; then
  echo "INFO this is $OS/$ARCH: the Linux installer is verified (it is what the GCP image uses) but not installed here."
  echo "INFO install the ORCA 6.1.0 build for this platform, then run check_toolchain.py --orca /ABS/PATH/orca"
  exit 0
fi

if [ -e "$PREFIX/orca" ]; then
  echo "INFO $PREFIX/orca already exists; skipping the installer and checking it"
else
  mkdir -p "$PREFIX"
  sh "$INSTALLER" -- -p "$PREFIX"
fi
[ -x "$PREFIX/orca" ] || { echo "FAIL $PREFIX/orca missing after install"; exit 1; }
BIN="$(sha256 "$PREFIX/orca")"
[ "$BIN" = "$WANT_BIN" ] || { echo "FAIL installed orca sha256 $BIN != locked $WANT_BIN"; exit 1; }
echo "PASS installed orca binary matches lock"

python3 - "$PREFIX" "$NAME" "$SHA" "$BIN" <<'PY'
import datetime, json, platform, sys
prefix, name, sha, binsha = sys.argv[1:]
json.dump({"orca_version": "6.1.0", "installer": name, "installer_sha256": sha,
           "orca_binary": prefix + "/orca", "orca_binary_sha256": binsha,
           "host": platform.node(), "installed_utc": datetime.datetime.utcnow().isoformat() + "Z"},
          open(prefix + "/chemflow_orca_receipt.json", "w"), indent=2)
PY
echo "DONE ORCA 6.1.0 at $PREFIX. ORCA needs Open MPI 4.1.8 on PATH for parallel runs:"
echo "     export PATH=$PREFIX:\$PATH; export LD_LIBRARY_PATH=$PREFIX:\${LD_LIBRARY_PATH:-}"
echo "Next: python3 tools/toolchain/check_toolchain.py --orca $PREFIX/orca"
