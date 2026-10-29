import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path(__file__).resolve().parent
RUN = BASE.parent

sim = pd.read_csv(BASE / "L3C_ATMOD7_WS_puro.csv")
exp = pd.read_csv(RUN / "L3C_datos_experimentales.csv")

assert len(sim) == len(exp)
assert np.allclose(
    sim["p_GeV"],
    exp["p_representativo_GeV"]
)

p = sim["p_GeV"].to_numpy(dtype=float)
I_sim = sim["Phi_total"].to_numpy(dtype=float)
sigma_sim = sim["sigma_total"].to_numpy(dtype=float)

I_exp = exp["Phi_m2"].to_numpy(dtype=float)
sigma_exp = exp["sigma_total_m2"].to_numpy(dtype=float)

validos = (
    np.isfinite(p) &
    np.isfinite(I_sim) &
    np.isfinite(I_exp) &
    np.isfinite(sigma_sim) &
    np.isfinite(sigma_exp) &
    (p > 0) &
    (I_sim > 0) &
    (I_exp > 0) &
    (sigma_sim >= 0) &
    (sigma_exp >= 0)
)

p = p[validos]
I_sim = I_sim[validos]
I_exp = I_exp[validos]
sigma_sim = sigma_sim[validos]
sigma_exp = sigma_exp[validos]

R = I_sim / I_exp

sigma_R = np.sqrt(
    (sigma_sim / I_exp)**2 +
    (I_sim * sigma_exp / I_exp**2)**2
)

# Estilo comun
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "legend.fontsize": 11,
    "figure.dpi": 130
})

# ==========================================
# FIGURA 1: ESPECTRO DIFERENCIAL
# ==========================================

fig, ax = plt.subplots(figsize=(9, 6))

ax.errorbar(
    p, I_exp,
    yerr=sigma_exp,
    fmt="o",
    color="black",
    markersize=5,
    capsize=2,
    label="L3+C"
)

ax.errorbar(
    p, I_sim,
    yerr=sigma_sim,
    fmt="s-",
    color="crimson",
    linewidth=1.8,
    markersize=5,
    capsize=3,
    label="CORSIKA-UrQMD (WS)"
)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"Momento del muon $p_\mu$ (GeV/c)")
ax.set_ylabel(
    r"Intensidad diferencial "
    r"[$m^{-2}\,s^{-1}\,sr^{-1}\,(GeV/c)^{-1}$]"
)

ax.set_title(
    "L3+C vs CORSIKA-UrQMD - ATMOD7"
)

ax.grid(True, which="both", alpha=0.22)
ax.legend()
fig.tight_layout()

fig.savefig(
    BASE / "L3C_WS_FINAL_espectro.png",
    dpi=300
)
fig.savefig(
    BASE / "L3C_WS_FINAL_espectro.pdf"
)
plt.close(fig)

# ==========================================
# FIGURA 2: COCIENTE
# ==========================================

fig, ax = plt.subplots(figsize=(9, 5.5))

ax.axhspan(
    0.85, 1.15,
    color="gray",
    alpha=0.16,
    label=r"Franja $\pm 15\%$"
)

ax.axhline(
    1.0,
    color="black",
    linewidth=1.2,
    label="Acuerdo exacto"
)

ax.errorbar(
    p, R,
    yerr=sigma_R,
    fmt="s",
    color="crimson",
    markersize=6,
    capsize=3,
    elinewidth=1.5,
    label="CORSIKA / L3+C"
)

ax.set_xscale("log")
ax.set_xlabel(r"Momento del muon $p_\mu$ (GeV/c)")
ax.set_ylabel(
    r"$I_\mu^{\mathrm{CORSIKA}}/"
    r"I_\mu^{\mathrm{L3+C}}$"
)

ax.set_title(
    "Cociente simulacion / experimento"
)

ax.grid(True, which="both", alpha=0.22)
ax.legend()
fig.tight_layout()

fig.savefig(
    BASE / "L3C_WS_FINAL_cociente.png",
    dpi=300
)
fig.savefig(
    BASE / "L3C_WS_FINAL_cociente.pdf"
)
plt.close(fig)

# Resultados numericos
tabla = pd.DataFrame({
    "p_GeV": p,
    "I_sim": I_sim,
    "sigma_sim": sigma_sim,
    "I_exp": I_exp,
    "sigma_exp": sigma_exp,
    "R": R,
    "sigma_R": sigma_R,
    "desviacion_pct": 100*(R-1)
})

tabla.to_csv(
    BASE / "L3C_WS_FINAL_resultados.csv",
    index=False
)

print("=== FIGURAS FINALES GENERADAS ===")
print("L3C_WS_FINAL_espectro.png / .pdf")
print("L3C_WS_FINAL_cociente.png / .pdf")
print("L3C_WS_FINAL_resultados.csv")
print("Puntos validos:", len(p))
print(
    "Dentro de +/-15 %:",
    np.sum(np.abs(R-1) <= 0.15)
)
print(
    "Compatibles con R=1 a 1 sigma:",
    np.sum(np.abs(R-1) <= sigma_R)
)
