#!/usr/bin/env python3
"""HydroSink progress dashboard: P1/P2 totals and one row per machine (this PC + Google Cloud VMs).

Machine rows are read live from each VM over SSH (8 s timeout); if a VM cannot be reached, its row falls back
to the copy synced every 10 min (analysis/logs/gcp/<vm>/) and is marked "(synced copy)".
'Results here' = marker tables present on this PC (cloud results arrive with the 10-min sync).
"""
import datetime as dt
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
A = os.path.join(ROOT, "analysis")
G = "/home/angsuman/extra_spac/software_libs/google-cloud-sdk/bin/gcloud"
PROJ, ZONE = os.environ.get("HS_GCP_PROJECT", ""), os.environ.get("HS_GCP_ZONE", "us-east4-c")
VMS = []  # all VMs deleted 2026-11-17 after final sync
KEY = os.path.expanduser("~/.ssh/google_compute_engine")
KH = os.path.expanduser("~/.ssh/google_compute_known_hosts")
VM_A = "/home/angsuman/extra_spac/Metagenomic Paper/analysis"
NOW = dt.datetime.now(dt.timezone.utc)


def lines(p):
    try:
        return [l.rstrip("\n") for l in open(p) if l.strip()]
    except OSError:
        return []


def done_have(runs):
    return sum(os.path.getsize(f) > 0 for r in runs
               if os.path.exists(f := os.path.join(A, "results/per_sample", f"{r}.markers.tsv")))


def parse_done(rows):
    """done.tsv rows -> list of (run, aware datetime); handles +05:30 and +00:00 stamps."""
    out = []
    for l in rows:
        f = l.split("\t")
        if len(f) >= 2:
            try:
                out.append((f[0], dt.datetime.fromisoformat(f[1])))
            except ValueError:
                pass
    return out


def rate(done, runs, hours=1.0):
    s = NOW - dt.timedelta(hours=hours)
    return sum(1 for r, t in done if r in runs and t >= s)


def vm_ips():
    try:
        o = subprocess.run([G, "compute", "instances", "list", "--project", PROJ,
                            "--format=value(name,networkInterfaces[0].accessConfigs[0].natIP,status)"],
                           capture_output=True, text=True, timeout=30).stdout
        return {f[0]: (f[1] if len(f) > 2 else "", f[-1]) for f in (l.split("\t") for l in o.splitlines()) if f}
    except Exception:
        return {}


def vm_live(name, ip):
    cmd = (f'D="{VM_A}"; echo "@@DONE"; cat "$D/logs/done.tsv" 2>/dev/null; echo "@@FAIL"; '
           f'cut -f1 "$D/logs/failures.tsv" 2>/dev/null | sort -u | wc -l; echo "@@DM"; pgrep -x diamond | wc -l; '
           f'echo "@@CPU"; top -bn2 -d1 | grep "Cpu(s)" | tail -1 | awk \'{{print 100-$8}}\'; '
           f'echo "@@MEM"; free -g | awk \'/Mem/{{print $3"/"$2}}\'; echo "@@END"; test -f "$D/logs/ALL_DONE" && echo yes')
    try:
        o = subprocess.run(["ssh", "-i", KEY, "-o", "StrictHostKeyChecking=no", "-o", f"UserKnownHostsFile={KH}",
                            "-o", "ConnectTimeout=8", "-o", "BatchMode=yes", f"angsuman@{ip}", cmd],
                           capture_output=True, text=True, timeout=25).stdout
    except Exception:
        return None
    if "@@END" not in o:
        return None
    sec, cur = {}, None
    for l in o.splitlines():
        if l.startswith("@@"):
            cur = l[2:]; sec[cur] = []
        elif cur:
            sec[cur].append(l)
    g = lambda k: (sec.get(k) or ["?"])[0].strip()
    return dict(done=parse_done(sec.get("DONE", [])), fail=g("FAIL"), dm=g("DM"), cpu=g("CPU"),
                mem=g("MEM") + " GB", finished=bool(sec.get("END")))


def fmt_eta(left, r):
    if left <= 0:
        return "finished"
    return f"{left / r:.1f} h" if r > 0 else "n/a"


def main():
    p2 = lines(os.path.join(A, "inventory/p2_queue.txt"))
    p1 = lines(os.path.join(A, "inventory/p1_queue.txt"))
    claims = set(lines(os.path.join(A, "inventory/p2_remote_claims.txt")))
    print(f"HydroSink progress   {dt.datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"P1 : {done_have(p1)}/{len(p1)}")
    p2have = done_have(p2)
    print(f"P2 : {p2have}/{len(p2)} results on this PC   ({len(p2) - p2have} to go)")
    print()
    hdr = f"{'machine':<14}{'done/assigned':>14}{'last hour':>11}{'ETA':>9}{'searches':>10}{'CPU':>6}{'RAM':>11}{'fail':>6}  note"
    print(hdr); print("-" * len(hdr))

    # this PC
    local_runs = [r for r in p2 if r not in claims]
    ldone = parse_done(lines(os.path.join(A, "logs/done.tsv")))
    lset = set(local_runs)
    lhave = done_have(local_runs)
    lr = rate(ldone, lset)
    dm = subprocess.run(["pgrep", "-xc", "diamond"], capture_output=True, text=True).stdout.strip() or "0"
    try:
        top = subprocess.run(["top", "-bn2", "-d1"], capture_output=True, text=True).stdout
        cpu = f"{100 - float([l for l in top.splitlines() if 'Cpu(s)' in l][-1].split(',')[3].split()[0]):.0f}"
    except Exception:
        cpu = "?"
    mem = subprocess.run(["free", "-g"], capture_output=True, text=True).stdout.splitlines()[1].split()
    act = subprocess.run(["systemctl", "--user", "is-active", "hs-p2"], capture_output=True, text=True).stdout.strip()
    fails = len({l.split("\t")[0] for l in lines(os.path.join(A, "logs/failures.tsv"))
                 if "offsets" in l and l.split("\t")[0] in lset})
    print(f"{'this PC':<14}{f'{lhave}/{len(local_runs)}':>14}{lr:>11}{fmt_eta(len(local_runs) - lhave, lr):>9}"
          f"{dm:>10}{cpu + '%':>6}{mem[2] + '/' + mem[1] + ' GB':>11}{fails:>6}  job {act}")
    tot_rate, etas = lr, [((len(local_runs) - lhave) / lr) if lr else (0 if lhave >= len(local_runs) else None)]

    ips = vm_ips()
    with ThreadPoolExecutor(max(1, len(VMS))) as ex:
        live = dict(zip(VMS, ex.map(lambda v: vm_live(v, ips.get(v, ("", ""))[0]) if ips.get(v, ("", ""))[0] else None, VMS)))
    for v in VMS:
        runs = lines(os.path.join(A, f"inventory/p2_vm_{v}.txt")); rs = set(runs)
        info, note = live[v], ""
        if info is None:
            d = os.path.join(A, "logs/gcp", v)
            info = dict(done=parse_done(lines(os.path.join(d, "done.tsv"))),
                        fail=str(len({l.split('\t')[0] for l in lines(os.path.join(d, 'failures.tsv'))})),
                        dm="-", cpu="-", mem="-", finished=os.path.exists(os.path.join(d, "ALL_DONE")))
            note = f"(synced copy; VM {ips.get(v, ('', 'unknown'))[1].lower() or 'unknown'})"
        dn = len({r for r, _ in info["done"] if r in rs})
        r1 = rate(info["done"], rs)
        tot_rate += r1
        etas.append(((len(runs) - dn) / r1) if r1 else (0 if dn >= len(runs) else None))
        if info["finished"]:
            note = (note + " all passes finished").strip()
        cpu = info["cpu"] if info["cpu"] in ("-", "?") else f"{float(info['cpu']):.0f}"
        print(f"{v:<14}{f'{dn}/{len(runs)}':>14}{r1:>11}{fmt_eta(len(runs) - dn, r1):>9}{info['dm']:>10}"
              f"{cpu + ('%' if cpu not in ('-', '?') else ''):>6}{info['mem']:>11}{info['fail']:>6}  {note}")
    print("-" * len(hdr))
    known = [e for e in etas if e is not None]
    eta = f"{max(known):.1f} h (slowest machine)" if known and len(known) == len(etas) else "n/a (some machine has no rate yet)"
    print(f"{'all machines':<14}{'':>14}{tot_rate:>11}   P2 finishes in about {eta}")
    sync = lines(os.path.join(A, "logs/gcp_sync.log"))
    print(f"\ncloud->PC copy: every 10 min; last: {sync[-1][:25] if sync else 'never'}"
          f"   (service hs-gcp-sync {subprocess.run(['systemctl', '--user', 'is-active', 'hs-gcp-sync'], capture_output=True, text=True).stdout.strip()})")


if __name__ == "__main__":
    main()
