#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if ! command -v "${CXX:-g++}" >/dev/null 2>&1 || ! command -v make >/dev/null 2>&1 ||
   { [[ ! -f /usr/include/NTL/lzz_pE.h && ! -f /usr/local/include/NTL/lzz_pE.h ]] && [[ -z ${NTL_PREFIX:-} ]]; }; then
  if [[ ${INSTALL_DEPS:-1} == 1 ]]; then bash install_deps.sh; else
    printf '%s\n' 'Missing build dependency. Run bash install_deps.sh or set NTL_PREFIX.' >&2
    exit 1
  fi
fi
make -j"${BUILD_JOBS:-2}" all
