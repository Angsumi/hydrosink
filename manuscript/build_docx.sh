#!/usr/bin/env bash
# Build the Word version of the manuscript (figures embedded after the legends). Needs pandoc >= 3 (PANDOC=...).
set -euo pipefail
cd "$(dirname "$0")/.."
# pandoc: $PANDOC, else on PATH, else the project's 'docs' conda env
DOCS_PANDOC=/home/angsuman/extra_spac/software_libs/conda_envs/docs/bin/pandoc
PANDOC=${PANDOC:-$(command -v pandoc || echo "$DOCS_PANDOC")}
[ -x "$PANDOC" ] || { echo "pandoc not found: set PANDOC=/path/to/pandoc"; exit 1; }
V=$(grep -oE 'Draft v[0-9.]+' manuscript/manuscript.md | head -n 1 | cut -c8-)
TITLE=$(sed -n '1s/^# //p' manuscript/manuscript.md)
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
tail -n +2 manuscript/manuscript.md > "$WORK/manuscript.md"
cat >> "$WORK/manuscript.md" <<'EOM'

## Figures
![Figure 1. Study design.](figures/fig1_study_design.png){width=16cm}

![Figure 2. Within-animal foregut-versus-hindgut contrast (P1).](figures/fig2_p1_within_animal.png){width=16cm}

![Figure 3. Faecal ALR_acet across vertebrates (P2).](figures/fig3_p2_cross_host.png){width=16cm}

![Figure 4. Proportionality between mcrA and acsB within studies (P3).](figures/fig4_p3_tradeoff.png){width=11cm}

![Figure S1. Marker detection by fermentation site.](figures/figS1_detection_by_site.png){width=16cm}

![Figure S2. PRISMA flow diagram.](figures/figS2_prisma_flow.png){width=13cm}
EOM
"$PANDOC" "$WORK/manuscript.md" -f markdown -t docx --metadata "title=$TITLE" \
  --resource-path=. -o "$WORK/manuscript.docx"
python3 - "$WORK/manuscript.docx" "manuscript/HydroSink_manuscript_v${V}.docx" <<'PYDOCX'
import sys
from xml.dom import minidom
from zipfile import ZipFile

word = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
with ZipFile(sys.argv[1]) as source, ZipFile(sys.argv[2], "w") as output:
    for item in source.infolist():
        data = source.read(item)
        if item.filename == "word/styles.xml":
            doc = minidom.parseString(data)
            for style in doc.getElementsByTagNameNS(word, "style"):
                style_id = style.getAttributeNS(word, "styleId")
                if style_id != "Title" and not style_id.startswith("Heading"):
                    continue
                for color in style.getElementsByTagNameNS(word, "color"):
                    for attr in ("themeColor", "themeTint", "themeShade"):
                        if color.hasAttributeNS(word, attr):
                            color.removeAttributeNS(word, attr)
                    color.setAttributeNS(word, "w:val", "000000")
                for border in style.getElementsByTagNameNS(word, "pBdr"):
                    border.parentNode.removeChild(border)
            data = doc.toxml(encoding="utf-8")
        output.writestr(item, data)
PYDOCX
echo "wrote manuscript/HydroSink_manuscript_v${V}.docx"
