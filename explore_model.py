from cobra.io import load_json_model
from pathlib import Path

model = load_json_model(str(Path("models/iJN678.json")))

# Print ALL reactions containing isoprenoid-related metabolites
# BiGG uses different naming -- search by metabolite, not reaction name
isoprenoid_metabolites = ["ipdp", "dmpp", "frdp", "grdp", "pyr", "g3p"]

print("Reactions involving isoprenoid metabolites:")
for met_id in isoprenoid_metabolites:
    try:
        met = model.metabolites.get_by_id(f"{met_id}_c")
        print(f"\n  {met_id} ({met.name}):")
        for rxn in met.reactions:
            flux = 0.0
            print(f"    {rxn.id:30s} {rxn.reaction}")
    except KeyError:
        print(f"\n  {met_id}: NOT FOUND in model")