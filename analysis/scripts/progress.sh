#!/usr/bin/env bash
# Progress dashboard for all machines (this PC + Google Cloud VMs). Live view: watch -n 60 bash analysis/scripts/progress.sh
exec python3 -I "$(cd "$(dirname "$0")" && pwd)/hs_dashboard.py"
