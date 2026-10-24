#!/usr/bin/env python3
"""Compile the completed literature workflow's JSON into review artifacts."""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = Path('/home/ubuntu/.manus-jobs/44faf84cb52e/output.txt')
raw = SOURCE.read_text()
# The addon output is one JSON object on its first line.
data = json.loads(raw.splitlines()[0])
papers = data['papers']
(ROOT / 'literature_verified.json').write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')

fields = ['id','citation','doi_or_identifier','primary_url','research_question','methods_and_model','main_findings','limitations','relationship_to_projects','full_text_status','evidence_stance']
with (ROOT / 'literature_comparison.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for p in papers:
        writer.writerow({k: p.get(k, '') for k in fields})

area = {
    '1':'Infinite-layer nickelate experiments',
    '2':'Infinite-layer nickelate experiments',
    '3':'Infinite-layer nickelate experiments',
    '4':'Bulk nickelate electronic structure',
    '5':'Bulk nickelate electronic structure',
    '6':'Bulk nickelate electronic structure',
    '7':'Oxide-interface transport',
    '8':'Oxide-interface electronic reconstruction',
    '9':'Oxide-interface electrostatics',
    '10':'Correlated-disorder critical theory',
    '11':'Correlated lattice percolation',
    '12':'Correlated lattice percolation',
}
lines = [
    '# Verified Literature Review', '',
    'This review summarizes 12 independently checked primary papers selected for the nickelate experiments, bulk/interface electronic structure, oxide-interface electrostatics, and correlated percolation. The two citations supplied in the user manuscript are reviewed separately in the main diagnosis as references [13] and [14].',
    '',
    '## Synthesis', '',
    data['synthesis']['literature_synthesis'], '',
    '## Paper-by-paper evidence', '',
]
for p in papers:
    lines += [
        f"### [{p['id']}] {p['citation']}",
        f"**Area:** {area.get(p['id'], 'Related background')}. **DOI/identifier:** {p.get('doi_or_identifier','')}.",
        f"**Research question:** {p.get('research_question','')}",
        f"**Methods/model:** {p.get('methods_and_model','')}",
        f"**Main findings:** {p.get('main_findings','')}",
        f"**Limitations and evidence stance:** {p.get('limitations','')} {p.get('evidence_stance','')}",
        f"**Relevance to the projects:** {p.get('relationship_to_projects','')}",
        f"**Access checked:** {p.get('full_text_status','')}",
        f"**Primary record:** {p.get('primary_url','')}",
        '',
    ]
lines += ['## Scope and access caveats', '']
lines.extend('- ' + x for x in data['synthesis'].get('source_caveats', []))
(ROOT / 'literature_review.md').write_text('\n'.join(lines) + '\n')
print(f"Compiled {len(papers)} verified papers; wrote literature_review.md, literature_comparison.csv, and literature_verified.json")
