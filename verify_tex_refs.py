"""Static LaTeX audit for cross-references and equation numbering.

The check is intentionally dependency-free. It verifies the invariants that are
most likely to break while editing the manuscript:
- no duplicate labels;
- no undefined cross-references;
- every eq:* label is actually referenced;
- every numbered display has a label;
- no hard-coded prose references such as "equations (7.20)--(7.21)";
- no manual \\tag numbering;
- balanced LaTeX environments at the begin/end-count level;
- all \\input files and bibliography citation keys exist;
- equation numbers are section-based.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
TEX_FILES = [ROOT / "main.tex", *sorted((ROOT / "chapter").glob("*.tex"))]

texts = {path: path.read_text(encoding="utf-8") for path in TEX_FILES}
joined = "\n".join(f"% FILE: {path.relative_to(ROOT)}\n{text}" for path, text in texts.items())

errors: list[str] = []

labels = re.findall(r"\\label\{([^}]+)\}", joined)
label_counts = Counter(labels)
for label, count in sorted(label_counts.items()):
    if count > 1:
        errors.append(f"duplicate label: {label} ({count} occurrences)")

refs = re.findall(r"\\(?:ref|eqref|cref|Cref)\{([^}]+)\}", joined)
for ref in sorted(set(refs)):
    if ref not in label_counts:
        errors.append(f"undefined reference: {ref}")

for label in sorted(set(labels)):
    if label.startswith("eq:") and label not in refs:
        errors.append(f"unused equation label: {label}")

numbered_envs = ("equation", "align", "gather", "multline")
for path, text in texts.items():
    rel = path.relative_to(ROOT)
    for env in numbered_envs:
        pattern = re.compile(
            rf"\\begin\{{{env}\}}([\s\S]*?)\\end\{{{env}\}}"
        )
        for match in pattern.finditer(text):
            if r"\label{" not in match.group(1):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: numbered {env} has no label")

hard_number = re.compile(
    r"\b(?:equations?|Eq\.?|Eqs\.?)\s*\(?\d+(?:\.\d+)+"
    r"(?:\s*--\s*\(?\d+(?:\.\d+)+\)?)?",
    re.IGNORECASE,
)
for path, text in texts.items():
    rel = path.relative_to(ROOT)
    for lineno, line in enumerate(text.splitlines(), start=1):
        if hard_number.search(line):
            errors.append(
                f"{rel}:{lineno}: hard-coded equation-number reference: {line.strip()}"
            )
        if r"\tag{" in line or r"\tag*{" in line:
            errors.append(f"{rel}:{lineno}: manual equation tag is not allowed")

begins = Counter(re.findall(r"\\begin\{([^}]+)\}", joined))
ends = Counter(re.findall(r"\\end\{([^}]+)\}", joined))
for env in sorted(set(begins) | set(ends)):
    if begins[env] != ends[env]:
        errors.append(
            f"environment count mismatch for {env}: "
            f"begin={begins[env]}, end={ends[env]}"
        )

main_text = texts[ROOT / "main.tex"]
if r"\numberwithin{equation}{section}" not in main_text:
    errors.append("main.tex must use section-based equation numbering")

for raw in re.findall(r"\\input\{([^}]+)\}", main_text):
    path = ROOT / (raw if raw.endswith(".tex") else raw + ".tex")
    if not path.exists():
        errors.append(f"missing input file: {path.relative_to(ROOT)}")

bib_path = ROOT / "references.bib"
bib_text = bib_path.read_text(encoding="utf-8")
bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib_text))
cite_chunks = re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}", joined)
cite_keys = {
    key.strip()
    for chunk in cite_chunks
    for key in chunk.split(",")
    if key.strip()
}
for key in sorted(cite_keys - bib_keys):
    errors.append(f"missing bibliography key: {key}")

if errors:
    print("LaTeX audit FAILED:")
    for item in errors:
        print(f"  - {item}")
    sys.exit(1)

print(
    "LaTeX audit passed: "
    f"{len(labels)} labels, {len(refs)} cross-references, "
    f"{len(cite_keys)} citation keys; numbering invariants are clean."
)
