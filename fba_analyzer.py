import cobra
from cobra.io import load_json_model
from cobra.flux_analysis import flux_variability_analysis

# 1. Load your PCC 11901 proxy model
model = load_json_model("models/iJN678.json")

# 2. Inspect the model -- what reactions exist?
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")

# 3. Baseline FBA -- maximize biomass (default objective)
with model:
    baseline = model.optimize()
    print(f"Baseline biomass flux: {baseline.objective_value:.4f}")

# 4. Find MEP pathway reactions -- your precursor supply chain
mep_keywords = ["dxs", "dxr", "ippi", "fpps", "mep", "dxp"]
mep_reactions = [
    r for r in model.reactions
    if any(k in r.id.lower() for k in mep_keywords)
]
print(f"\nMEP pathway reactions found: {len(mep_reactions)}")
for r in mep_reactions:
    print(f"  {r.id:25s} flux={baseline.fluxes.get(r.id, 0):.4f}")

# 5. Add SSL-3 reaction (botryococcene synthase)
#    FPP + FPP → botryococcene  (simplified -- SSL-1 makes PSPP, SSL-3 makes botryococcene)
#    For FBA purposes, model as: 2 FPP → botryococcene
with model:
    # Add botryococcene as a metabolite
    botryococcene = cobra.Metabolite(
        id="botryococcene_c",
        name="Botryococcene",
        compartment="c",
        formula="C30H50"
    )

    # Add SSL-3 reaction
    ssl3 = cobra.Reaction("SSL3")
    ssl3.name = "SSL-3: Botryococcene synthase"
    ssl3.lower_bound = 0
    ssl3.upper_bound = 1000

    # Find FPP metabolite in the model
    fpp = model.metabolites.get_by_id("frdp_c")  # farnesyl-PP in BiGG notation
    ssl3.add_metabolites({
        fpp: -2,                  # consumes 2 FPP
        botryococcene: 1          # produces botryococcene
    })
    model.add_reactions([ssl3])

    # Add exchange reaction so botryococcene can leave the system
    exchange = cobra.Reaction("EX_botryococcene")
    exchange.add_metabolites({botryococcene: -1})
    exchange.lower_bound = 0
    exchange.upper_bound = 1000
    model.add_reactions([exchange])

    # 6. Maximize botryococcene production
    model.objective = "EX_botryococcene"
    sol = model.optimize()
    print(f"\nMax botryococcene flux: {sol.objective_value:.4f} mmol/gDW/h")
    print(f"Status: {sol.status}")

    # 7. FVA -- find bottlenecks
    fva = flux_variability_analysis(
        model,
        reaction_list=mep_reactions,
        fraction_of_optimum=0.9
    )
    print("\nMEP pathway flux variability (at 90% optimum):")
    print(fva)