#!/usr/bin/env bash
# One-time setup of a Google Cloud worker VM for P2 (amendment v1.2 protocol).
# Recreates the local directory layout so 02_fetch_quantify.sh runs unchanged, with the same tool versions:
# diamond 2.1.10, fastp 0.23.4, python 3.11, sra-tools 3.1.1 (bioconda / conda-forge).
set -euo pipefail
SW=/home/angsuman/extra_spac/software_libs
mkdir -p "$SW/bin" "$SW/conda_envs"
if [ ! -x "$SW/bin/micromamba" ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$SW" bin/micromamba
fi
export MAMBA_ROOT_PREFIX=$SW/mamba_root
MM="$SW/bin/micromamba"
[ -x "$SW/conda_envs/hydrosink/bin/diamond" ] || "$MM" create -y -q -p "$SW/conda_envs/hydrosink" -c conda-forge -c bioconda \
  diamond=2.1.10 fastp=0.23.4 python=3.11 pigz
[ -x "$SW/conda_envs/sratools/bin/vdb-dump" ] || "$MM" create -y -q -p "$SW/conda_envs/sratools" -c conda-forge -c bioconda sra-tools=3.1.1
"$SW/conda_envs/hydrosink/bin/diamond" version
"$SW/conda_envs/hydrosink/bin/fastp" --version
"$SW/conda_envs/hydrosink/bin/python3.11" --version
"$SW/conda_envs/sratools/bin/vdb-dump" --version | grep -m1 version
