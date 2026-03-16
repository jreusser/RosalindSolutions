"""
Add SSL-1/2/3 (B. braunii botryococcene pathway) to iJN678
and run FBA + FVA to identify flux bottlenecks.

Key finding from explore_model.py:
- FPP = frdp_c in BiGG notation
- Native competitor: SQLS (2x frdp_c → squalene)
- Supply chain: DXPS → DXPRIi → ... → IPDPS_syn → IPDDI → DMATT → GRTT → frdp_c
"""
import cobra
from cobra.io import load_json_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import production_envelope
from pathlib import Path
import numpy as np
import pandas as pd

MODEL_PATH = Path("models/iJN678.json")

# BiGG metabolite IDs confirmed from explore_model.py
FPP_ID    = "frdp_c"
NADPH_ID  = "nadph_c"
NADP_ID   = "nadp_c"
PPI_ID    = "ppi_c"
H_ID      = "h_c"

def load_model() -> cobra.Model:
    model = load_json_model(str(MODEL_PATH))
    print(f"Model: {model.id}")
    print(f"Reactions:   {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    return model

def add_ssl_pathway(model: cobra.Model) -> cobra.Model:
    """
    Add SSL-1/2/3 reactions to the model.
    
    SSL-1: FPP + FPP → PSPP + PPi        (presqualene diphosphate synthase)
    SSL-2: PSPP + NADPH → squalene + PPi  (squalene synthase — membrane)
    SSL-3: PSPP + NADPH → botryococcene   (botryococcene synthase — product)
    
    Stoichiometry from: Niehaus et al. 2011, PNAS
    """

    # --- Metabolites ---
    pspp = cobra.Metabolite(
        id="pspp_c",
        name="Presqualene diphosphate",
        compartment="c",
        formula="C30H52O7P2"
    )
    botryococcene = cobra.Metabolite(
        id="botryococcene_c",
        name="Botryococcene",
        compartment="c",
        formula="C30H50"
    )

    fpp   = model.metabolites.get_by_id(FPP_ID)
    nadph = model.metabolites.get_by_id(NADPH_ID)
    nadp  = model.metabolites.get_by_id(NADP_ID)
    ppi   = model.metabolites.get_by_id(PPI_ID)
    h     = model.metabolites.get_by_id(H_ID)

    # --- SSL-1: 2 FPP → PSPP + PPi ---
    ssl1 = cobra.Reaction("SSL1")
    ssl1.name       = "SSL-1 Presqualene diphosphate synthase (B. braunii)"
    ssl1.lower_bound = 0
    ssl1.upper_bound = 1000
    ssl1.add_metabolites({
        fpp:  -2,
        pspp:  1,
        ppi:   1,
    })

    # --- SSL-2: PSPP + NADPH + H → squalene + NADP + PPi ---
    # squalene (sql_c) already exists in model via SQLS reaction
    sql = model.metabolites.get_by_id("sql_c")
    ssl2 = cobra.Reaction("SSL2")
    ssl2.name        = "SSL-2 Squalene synthase (B. braunii) — membrane hopanoids"
    ssl2.lower_bound = 0
    ssl2.upper_bound = 1000
    ssl2.add_metabolites({
        pspp:  -1,
        nadph: -1,
        h:     -1,
        sql:    1,
        nadp:   1,
        ppi:    1,
    })

    # --- SSL-3: PSPP + NADPH + H → botryococcene + NADP + PPi ---
    ssl3 = cobra.Reaction("SSL3")
    ssl3.name        = "SSL-3 Botryococcene synthase (B. braunii) — fuel product"
    ssl3.lower_bound = 0
    ssl3.upper_bound = 1000
    ssl3.add_metabolites({
        pspp:          -1,
        nadph:         -1,
        h:             -1,
        botryococcene:  1,
        nadp:           1,
        ppi:            1,
    })

    # --- Exchange reaction: botryococcene leaves the cell ---
    ex_bot = cobra.Reaction("EX_botryococcene_c")
    ex_bot.name        = "Botryococcene export"
    ex_bot.lower_bound = 0
    ex_bot.upper_bound = 1000
    ex_bot.add_metabolites({botryococcene: -1})

    model.add_reactions([ssl1, ssl2, ssl3, ex_bot])
    print(f"\nSSL-1/2/3 added. New reaction count: {len(model.reactions)}")
    return model

def run_analysis(model: cobra.Model) -> None:

    # --- Baseline: maximize biomass (unchanged objective) ---
    with model:
        baseline = model.optimize()
        native_sqls = baseline.fluxes.get("SQLS", 0.0)
        print(f"\n--- Baseline (maximize biomass) ---")
        print(f"Biomass flux:           {baseline.objective_value:.6f}")
        print(f"Native SQLS flux:       {native_sqls:.6f}  (squalene — membrane)")
        print(f"SSL1 flux:              {baseline.fluxes.get('SSL1', 0.0):.6f}")
        print(f"SSL3 flux:              {baseline.fluxes.get('SSL3', 0.0):.6f}")
        print(f"Botryococcene export:   {baseline.fluxes.get('EX_botryococcene_c', 0.0):.6f}")

    # --- Maximize botryococcene ---
    with model:
        model.objective = "EX_botryococcene_c"
        sol = model.optimize()
        print(f"\n--- Maximize botryococcene ---")
        print(f"Status:                 {sol.status}")
        print(f"Max botryococcene flux: {sol.objective_value:.6f} mmol/gDW/h")
        print(f"Biomass flux:           {sol.fluxes.get('BIOMASS_Cyanosys', 0.0):.6f}")
        print(f"Native SQLS flux:       {sol.fluxes.get('SQLS', 0.0):.6f}")
        print(f"SSL1 flux:              {sol.fluxes.get('SSL1', 0.0):.6f}")
        print(f"SSL2 flux:              {sol.fluxes.get('SSL2', 0.0):.6f}")
        print(f"SSL3 flux:              {sol.fluxes.get('SSL3', 0.0):.6f}")
        print(f"DXPS flux:              {sol.fluxes.get('DXPS', 0.0):.6f}")
        print(f"GRTT flux (→FPP):       {sol.fluxes.get('GRTT', 0.0):.6f}")

        # --- FVA on MEP + SSL reactions ---
        mep_ssl_reactions = [
            "DXPS", "DXPRIi", "MEPCT",
            "IPDDI", "DMATT", "GRTT",
            "SSL1", "SSL2", "SSL3",
            "SQLS",
        ]
        # filter to only reactions that exist in model
        valid = [r for r in mep_ssl_reactions if r in [x.id for x in model.reactions]]

        print(f"\n--- FVA at 90% optimum (botryococcene objective) ---")
        fva = flux_variability_analysis(
            model,
            reaction_list=valid,
            fraction_of_optimum=0.9
        )
        print(fva.to_string())

        # --- The key ratio: SSL2 vs SSL3 competition for PSPP ---
        ssl2_flux = sol.fluxes.get("SSL2", 0.0)
        ssl3_flux = sol.fluxes.get("SSL3", 0.0)
        total     = ssl2_flux + ssl3_flux
        if total > 0:
            print(f"\n--- PSPP Flux Split ---")
            print(f"SSL2 (→ squalene):      {ssl2_flux:.6f} ({ssl2_flux/total*100:.1f}%)")
            print(f"SSL3 (→ botryococcene): {ssl3_flux:.6f} ({ssl3_flux/total*100:.1f}%)")
            print(f"Total PSPP consumed:    {total:.6f}")


def production_envelope_analysis(model: cobra.Model) -> None:
    """
    Compute the biomass–botryococcene Pareto front.
    Sweeps minimum biomass constraint and finds max botryococcene at each point.
    """
    print(f"\n{'='*60}")
    print(f"PRODUCTION ENVELOPE: Biomass vs Botryococcene")
    print(f"{'='*60}")

    # Find biomass reaction ID
    biomass_rxn = [r for r in model.reactions if "BIOMASS" in r.id.upper()]
    if not biomass_rxn:
        print("ERROR: No biomass reaction found.")
        return
    biomass_id = biomass_rxn[0].id

    with model:
        # First get max biomass
        model.objective = biomass_id
        max_biomass = model.optimize().objective_value

        # Now sweep biomass from 0% to 100% and maximize botryococcene
        model.objective = "EX_botryococcene_c"
        n_points = 20
        fractions = np.linspace(0, 1, n_points + 1)
        results = []

        biomass_rxn_obj = model.reactions.get_by_id(biomass_id)

        for frac in fractions:
            min_biomass = frac * max_biomass
            biomass_rxn_obj.lower_bound = min_biomass
            sol = model.optimize()
            if sol.status == "optimal":
                bot_flux = sol.objective_value
                dxps_flux = sol.fluxes.get("DXPS", 0.0)
                ssl1_flux = sol.fluxes.get("SSL1", 0.0)
            else:
                bot_flux = 0.0
                dxps_flux = 0.0
                ssl1_flux = 0.0
            results.append({
                "biomass_min": min_biomass,
                "biomass_frac": frac,
                "botryococcene": bot_flux,
                "DXPS": dxps_flux,
                "SSL1": ssl1_flux,
            })

    df = pd.DataFrame(results)
    print(df.to_string(index=False, float_format="{:.6f}".format))

    # Highlight practical operating points
    print(f"\n--- Practical Operating Points ---")
    print(f"{'Biomass %':>10} {'Biomass flux':>14} {'Botryo flux':>14} {'% of max botryo':>16}")
    max_bot = df["botryococcene"].max()
    for target_pct in [0, 10, 25, 50, 75, 90, 100]:
        row = df.iloc[(df["biomass_frac"] - target_pct/100).abs().argsort().iloc[0]]
        pct_of_max = (row["botryococcene"] / max_bot * 100) if max_bot > 0 else 0
        print(f"{target_pct:>9}% {row['biomass_min']:>14.6f} {row['botryococcene']:>14.6f} {pct_of_max:>15.1f}%")

    return df


def carbon_yield_analysis(model: cobra.Model) -> None:
    """
    Estimate theoretical carbon yield of botryococcene.
    Botryococcene = C30H50 → 30 carbons per molecule.
    """
    print(f"\n{'='*60}")
    print(f"CARBON YIELD ANALYSIS")
    print(f"{'='*60}")

    with model:
        model.objective = "EX_botryococcene_c"
        sol = model.optimize()
        max_bot = sol.objective_value  # mmol/gDW/h

        # Find CO2 exchange (carbon source for photoautotroph)
        co2_exchanges = [r for r in model.reactions if "co2" in r.id.lower() and "EX_" in r.id]
        print(f"\nCO2 exchange reactions:")
        for r in co2_exchanges:
            flux = sol.fluxes.get(r.id, 0.0)
            print(f"  {r.id:30s} flux={flux:.6f}  bounds=({r.lower_bound}, {r.upper_bound})")

        # Photon exchange
        photon_exchanges = [r for r in model.reactions if "photon" in r.id.lower() or "hnu" in r.id.lower()]
        print(f"\nPhoton exchange reactions:")
        for r in photon_exchanges:
            flux = sol.fluxes.get(r.id, 0.0)
            print(f"  {r.id:30s} flux={flux:.6f}  bounds=({r.lower_bound}, {r.upper_bound})")

        # Carbon yield: C in product / C consumed
        # Botryococcene = C30H50 → 30 C per mol
        c_in_product = max_bot * 30  # mmol C / gDW / h
        print(f"\nCarbon in botryococcene: {max_bot:.6f} mmol/gDW/h × 30 C = {c_in_product:.4f} mmol-C/gDW/h")

        # Molecular weight of C30H50 = 30(12) + 50(1) = 410 g/mol
        mass_flux = max_bot * 410.0 / 1000  # g/gDW/h (convert mmol to mol * g/mol)
        print(f"Mass flux: {max_bot:.6f} mmol/gDW/h × 410 g/mol = {mass_flux:.6f} g/gDW/h")
        print(f"           = {mass_flux*1000:.4f} mg/gDW/h")


def engineering_targets(model: cobra.Model) -> None:
    """
    Identify overexpression, knockout, and cofactor engineering targets
    based on FBA shadow prices and reduced costs.
    """
    print(f"\n{'='*60}")
    print(f"ENGINEERING TARGETS")
    print(f"{'='*60}")

    with model:
        model.objective = "EX_botryococcene_c"
        sol = model.optimize()

        # Shadow prices: metabolites whose limited supply constrains botryococcene
        print(f"\n--- Shadow Prices (top constraints on botryococcene) ---")
        print(f"Positive = increasing this metabolite's availability increases objective")
        shadows = sol.shadow_prices.sort_values(ascending=False)
        # Show top 15 positive and top 15 negative
        top_pos = shadows[shadows > 1e-6].head(15)
        top_neg = shadows[shadows < -1e-6].tail(15)
        if len(top_pos) > 0:
            print(f"\nMost beneficial to increase supply:")
            for met_id, price in top_pos.items():
                met = model.metabolites.get_by_id(met_id)
                print(f"  {met_id:25s} {price:>10.6f}  ({met.name})")
        if len(top_neg) > 0:
            print(f"\nMost beneficial to increase removal:")
            for met_id, price in top_neg.items():
                met = model.metabolites.get_by_id(met_id)
                print(f"  {met_id:25s} {price:>10.6f}  ({met.name})")

        # Reduced costs: reactions that would benefit from altered bounds
        print(f"\n--- Reduced Costs (reaction bound constraints) ---")
        reduced = sol.reduced_costs
        # Reactions at lower bound with negative reduced cost → useful if knocked in/up
        promising_up = reduced[reduced < -1e-6].sort_values().head(15)
        if len(promising_up) > 0:
            print(f"\nOverexpression targets (at lower bound, negative reduced cost):")
            for rxn_id, cost in promising_up.items():
                rxn = model.reactions.get_by_id(rxn_id)
                flux = sol.fluxes.get(rxn_id, 0.0)
                print(f"  {rxn_id:25s} cost={cost:>10.6f}  flux={flux:.6f}  ({rxn.name[:50]})")

        # Single-gene knockout screening for the MEP+SSL pathway
        print(f"\n--- Single Reaction Knockout Impact ---")
        targets = [
            "DXPS", "DXPRIi", "MEPCT", "IPDDI", "DMATT", "GRTT",
            "SSL1", "SSL2", "SSL3", "SQLS",
        ]
        valid_targets = [r for r in targets if r in [x.id for x in model.reactions]]
        for rxn_id in valid_targets:
            with model:
                model.reactions.get_by_id(rxn_id).knock_out()
                ko_sol = model.optimize()
                if ko_sol.status == "optimal":
                    print(f"  KO {rxn_id:12s} → botryococcene = {ko_sol.objective_value:.6f} "
                          f"({ko_sol.objective_value/sol.objective_value*100:.1f}% of WT)")
                else:
                    print(f"  KO {rxn_id:12s} → {ko_sol.status}")

        # SQLS knockout + biomass-constrained production
        print(f"\n--- SQLS Knockout + Growth-Coupled Production ---")
        with model:
            model.reactions.get_by_id("SQLS").knock_out()
            # First check if biomass is still feasible
            model.objective = model.reactions.get_by_id(
                [r.id for r in model.reactions if "BIOMASS" in r.id.upper()][0]
            )
            bio_sol = model.optimize()
            print(f"  ΔSQLS max biomass: {bio_sol.objective_value:.6f} "
                  f"(vs WT: {0.063150:.6f})")
            # Now maximize botryococcene at 50% max biomass
            model.objective = "EX_botryococcene_c"
            biomass_rxn = [r for r in model.reactions if "BIOMASS" in r.id.upper()][0]
            biomass_rxn.lower_bound = 0.5 * bio_sol.objective_value
            coupled_sol = model.optimize()
            if coupled_sol.status == "optimal":
                print(f"  ΔSQLS + 50% biomass: botryococcene = {coupled_sol.objective_value:.6f} mmol/gDW/h")
            else:
                print(f"  ΔSQLS + 50% biomass: {coupled_sol.status}")


if __name__ == "__main__":
    model = load_model()
    model = add_ssl_pathway(model)
    run_analysis(model)
    production_envelope_analysis(model)
    carbon_yield_analysis(model)
    engineering_targets(model)