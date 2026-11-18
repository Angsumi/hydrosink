#!/usr/bin/env bash
# Build the Word version of the manuscript (figures embedded after the legends). Needs pandoc >= 3 (PANDOC=...).
set -euo pipefail
cd "$(dirname "$0")/.."
PANDOC=${PANDOC:-/home/angsuman/extra_spac/software_libs/conda_envs/docs/bin/pandoc}
V=$(grep -oE 'Draft v[0-9.]+' manuscript/manuscript.md | head -n 1 | cut -c8-)
T=$(mktemp --suffix=.md); trap 'rm -f "$T"' EXIT
cp manuscript/manuscript.md "$T"
cat >> "$T" <<'EOM'

## Figures
![Figure 1. Study design.](figures/fig1_study_design.png){width=16cm}

![Figure 2. Within-animal foregut-versus-hindgut contrast (P1).](figures/fig2_p1_within_animal.png){width=16cm}

![Figure 3. Faecal ALR_acet across vertebrates (P2).](figures/fig3_p2_cross_host.png){width=16cm}

![Figure 4. Proportionality between mcrA and acsB within studies (P3).](figures/fig4_p3_tradeoff.png){width=11cm}

![Figure S1. Marker detection by fermentation site.](figures/figS1_detection_by_site.png){width=16cm}

![Figure S2. PRISMA flow diagram.](figures/figS2_prisma_flow.png){width=13cm}
EOM
"$PANDOC" "$T" -f markdown -t docx --resource-path=. -o "manuscript/HydroSink_manuscript_v${V}.docx"
echo "wrote manuscript/HydroSink_manuscript_v${V}.docx"
