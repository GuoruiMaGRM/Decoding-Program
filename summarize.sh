#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
[[ $# -ge 1 ]] || { printf 'Usage: bash summarize.sh RUN_DIR [EXPECTED_SEEDS]\n' >&2; exit 2; }
bash build.sh
bin/hhy_bench --task aggregate --root "$1" --expected-seeds "${2:-20}" --bootstrap "${BOOTSTRAP:-2000}"
