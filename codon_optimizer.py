import tellurium as te

model_str = """
# PCC 11901 — SSL Botryococcene Pathway
# Minimal SBML model for FBA / dynamic simulation

var IPP, DMAPP, FPP, PSPP, squalene, botryococcene

# MEP pathway flux (simplified input)
J_MEP: -> IPP; k_mep * light

# IDI: IPP ⇌ DMAPP
J_IDI: IPP -> DMAPP; k_idi * IPP - k_idi_r * DMAPP

# FPPS: IPP + DMAPP -> FPP
J_FPPS: IPP + DMAPP -> FPP; k_fpps * IPP * DMAPP

# SSL-1: FPP + FPP -> PSPP
J_SSL1: 2 FPP -> PSPP; k_ssl1 * FPP^2

# SSL-2: PSPP -> squalene (essential — membrane hopanoids)
J_SSL2: PSPP -> squalene; k_ssl2 * PSPP

# SSL-3: PSPP -> botryococcene (PRODUCT)
J_SSL3: PSPP -> botryococcene; k_ssl3 * PSPP

# Rate constants (placeholder — replace with literature Km/Vmax)
k_mep   = 0.5
k_idi   = 1.0; k_idi_r = 0.1
k_fpps  = 0.8
k_ssl1  = 0.6
k_ssl2  = 0.3  # tuned low — just enough squalene for membranes
k_ssl3  = 0.7  # tuned high — maximize botryococcene

light = 1.0    # normalized photon flux
"""

r = te.loadAntimonyModel(model_str)

# Export to SBML
sbml = r.getSBML()
with open("pcc11901_ssl_pathway.xml", "w") as f:
    f.write(sbml)

# Simulate 24h
result = r.simulate(0, 24, 1000)
r.plot()