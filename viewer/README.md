# Remission replay viewer: a run's lab work and tissue next to the reference scientist's, same patient, same clock
Rebuild data (stdlib Python, about 50 KB per run, one `LABEL=JOB_DIR` per row): `python3 scripts/build_viewer.py "Reference scientist=/bulk2/runs/cellworld/jobs/e2e1-oracle" "GPT-6.1 Sol=/bulk2/runs/cellworld/jobs/ag1-sol-agent" [--mistakes results/run_analysis.jsonl]`
Open: `python3 -m http.server -d viewer 8377`, then http://localhost:8377/ (opening `viewer/index.html` from disk also works; nothing is fetched from the network).
Deep link: `#patient=dev-01&left=gpt-6-1-sol&right=reference-scientist&day=96`; `#view=overview` is the patients × runs table. Keys: space, ← →, shift+← →, Home, End.
Without `--mistakes`, the "Where it went wrong" callouts are derived from the episode itself (which clone regrew, what was sequenced, tested, left unused).
