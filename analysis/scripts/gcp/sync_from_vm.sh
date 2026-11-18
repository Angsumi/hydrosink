#!/usr/bin/env bash
# Copy finished P2 results from the Google Cloud worker VMs to this machine, every 10 minutes.
# Only runs listed in a VM's done.tsv are copied (their files are complete); existing local files are never
# overwritten. VM logs are mirrored to analysis/logs/gcp/<vm>/. Usage: sync_from_vm.sh [once]
set -uo pipefail
ROOT="/home/angsuman/extra_spac/Metagenomic Paper"
G=/home/angsuman/extra_spac/software_libs/google-cloud-sdk/bin/gcloud
ZONE=${HS_GCP_ZONE:-us-east4-c}; PROJ=${HS_GCP_PROJECT:?set HS_GCP_PROJECT}; VMS=${HS_VMS:-"hs-p2-c56"}
KEY=$HOME/.ssh/google_compute_engine
SSH="ssh -i $KEY -o StrictHostKeyChecking=no -o UserKnownHostsFile=$HOME/.ssh/google_compute_known_hosts -o ConnectTimeout=20"
while :; do
  for VM in $VMS; do
    D="$ROOT/analysis/logs/gcp/$VM"; mkdir -p "$D"
    IP=$("$G" compute instances describe "$VM" --zone $ZONE --project $PROJ --format='value(networkInterfaces[0].accessConfigs[0].natIP)' 2>/dev/null)
    [ -z "$IP" ] && { echo "$(date -Is) $VM not reachable"; continue; }
    R="angsuman@$IP:/home/angsuman/extra_spac/Metagenomic Paper/analysis"  # rsync 3.2+ protects args: no escaping
    rsync -a -e "$SSH" "$R/logs/" "$D/" 2>/dev/null
    if [ -f "$D/done.tsv" ]; then
      L=$(mktemp)
      cut -f1 "$D/done.tsv" | sort -u | while read -r r; do
        printf 'results/per_sample/%s.markers.tsv\nresults/per_sample/%s.fastp.json\nresults/hits/%s.m8.gz\n' "$r" "$r" "$r"
      done > "$L"
      rsync -a --ignore-existing --files-from="$L" -e "$SSH" "$R/" "$ROOT/analysis/" 2>/dev/null
      rm -f "$L"
    fi
    echo "$(date -Is) $VM synced; done there: $(cut -f1 "$D/done.tsv" 2>/dev/null | sort -u | wc -l)"
  done
  [ "${1:-}" = once ] && break
  sleep 600
done
