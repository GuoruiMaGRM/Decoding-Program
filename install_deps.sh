#!/usr/bin/env bash
set -euo pipefail
# Invoked only on the user's machine. No dependency installation was run by
# the assistant. Existing NTL installations can use NTL_PREFIX instead.
if (( EUID == 0 )); then elevate=(); else elevate=(sudo); fi
if command -v apt-get >/dev/null 2>&1; then
  "${elevate[@]}" apt-get update
  "${elevate[@]}" apt-get install -y g++ make libntl-dev libgmp-dev coreutils util-linux
elif command -v dnf >/dev/null 2>&1; then
  "${elevate[@]}" dnf install -y gcc-c++ make ntl-devel gmp-devel coreutils util-linux
elif command -v pacman >/dev/null 2>&1; then
  "${elevate[@]}" pacman -S --needed --noconfirm gcc make ntl gmp coreutils util-linux
else
  printf '%s\n' 'Install a C++17 compiler, GNU make, NTL, GMP and coreutils, or set NTL_PREFIX.' >&2
  exit 1
fi
