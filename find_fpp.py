from cobra.io import load_json_model
from pathlib import Path

model = load_json_model(str(Path("models/iJN678.json")))

# Search for FPP by name, not ID
hits = [
    m for m in model.metabolites
    if any(k in m.name.lower() for k in [
        "farnesyl", "fpp", "farnesyl-pp",
        "farnesyl diphosphate", "farnesyl pyrophosphate"
    ])
]

print(f"FPP candidates ({len(hits)} found):")
for m in hits:
    print(f"  ID: {m.id:20s} Name: {m.name}")
    for r in m.reactions:
        print(f"    → {r.id:25s} {r.reaction[:60]}")