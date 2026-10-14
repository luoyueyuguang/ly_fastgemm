#!/usr/bin/env bash
# Fetches the comparison baselines used by bench_gemm.cpp into third_party/.
#
# Eigen, OpenBLAS and libxsmm are unpacked from their distribution packages
# rather than installed, so nothing outside this directory is touched and the
# library itself keeps having no dependency at all.  Needs apt-get; run once.
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p third_party/debs third_party/root third_party/lib

packages=(libeigen3-dev libopenblas-dev libopenblas0-pthread libxsmm-dev)
for pkg in "${packages[@]}"; do
  (cd third_party/debs && apt-get download "$pkg")
done

for deb in third_party/debs/*.deb; do
  dpkg-deb -x "$deb" third_party/root
done

# link names the compiler expects
openblas=$(find third_party/root -name 'libopenblasp-r*.so' | head -1)
if [ -z "$openblas" ]; then
  echo "no OpenBLAS shared object found" >&2
  exit 1
fi
ln -sf "../${openblas#third_party/}" third_party/lib/libopenblas.so
ln -sf ../root/usr/lib/libxsmm.a third_party/lib/libxsmm.a

echo "third_party/ is ready; now run: make bench"
