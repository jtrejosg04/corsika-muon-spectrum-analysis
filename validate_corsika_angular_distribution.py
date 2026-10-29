from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

BASE = Path(__file__).resolve().parent
SALIDA = BASE / "resultados_angulares_conjunto"
SALIDA.mkdir(exist_ok=True)

N_THETA, N_PHI = 6, 12
P_MIN, P_MAX = 1.0, 5.0
THETA_MAX = 30.0

primarios = []
muones = []

for run in range(401, 422):
    p = pd.read_csv(BASE / f"catalogos/primarios_{run}.csv")
    m = pd.read_csv(BASE / f"catalogos/muones_{run}.csv")
    p["run_id"] = run
    m["run_id"] = run
    primarios.append(p)
    muones.append(m)

prim = pd.concat(primarios, ignore_index=True)
mu = pd.concat(muones, ignore_index=True)

N = len(prim)
assert N == 105000

claves = ["run_id", "shower_id"]
assert not prim.duplicated(claves).any()

mu = mu[
    (mu.nivel == 1) &
    (mu.p >= P_MIN) &
    (mu.p <= P_MAX) &
    (mu.theta_deg >= 0) &
    (mu.theta_deg < THETA_MAX) &
    (mu.cos_theta > 0)
].copy()

mu["phi_360"] = np.mod(mu.phi_deg, 360.0)
mu["peso"] = 1.0 / mu.cos_theta

indice_prim = pd.MultiIndex.from_frame(prim[claves])
indice_mu = pd.MultiIndex.from_frame(mu[claves])
ids = indice_prim.get_indexer(indice_mu)

assert np.all(ids >= 0)

theta_edges = np.linspace(0, THETA_MAX, N_THETA + 1)
phi_edges = np.linspace(0, 360, N_PHI + 1)

def distribucion(variable, edges):
    nb = len(edges) - 1
    indices = np.searchsorted(edges, variable, side="right") - 1
    indices = np.clip(indices, 0, nb - 1)

    # Acumulacion dispersa por cascada e intervalo.
    # Evita crear una matriz enorme de N x nb.
    claves_bin = ids * nb + indices

    suma = np.bincount(
        indices,
        weights=mu.peso.to_numpy(),
        minlength=nb
    )

    # Suma de cuadrados de la contribucion por cascada.
    orden = np.argsort(claves_bin)
    claves_ordenadas = claves_bin[orden]
    pesos_ordenados = mu.peso.to_numpy()[orden]

    _, inicios = np.unique(
        claves_ordenadas, return_index=True
    )
    aportes = np.add.reduceat(pesos_ordenados, inicios)
    bins_unicos = claves_ordenadas[inicios] % nb

    suma_cuadrados = np.bincount(
        bins_unicos,
        weights=aportes**2,
        minlength=nb
    )

    varianza = (
        N / (N - 1)
        * (suma_cuadrados - suma**2 / N)
    )
    varianza = np.maximum(varianza, 0)

    anchura = np.diff(edges)
    return suma / anchura, np.sqrt(varianza) / anchura

def comparar(y, sigma, referencia):
    validos = np.isfinite(y) & (sigma > 0)
    yy, ss, rr = y[validos], sigma[validos], referencia[validos]

    a = np.sum(yy * rr / ss**2) / np.sum(rr**2 / ss**2)
    chi = np.sum(((yy - a * rr) / ss)**2)
    gl = len(yy) - 1

    return a, chi, gl, chi2.sf(chi, gl)

yth, eth = distribucion(
    mu.theta_deg.to_numpy(), theta_edges
)
yph, eph = distribucion(
    mu.phi_360.to_numpy(), phi_edges
)

t = np.deg2rad(theta_edges)
dth = np.diff(theta_edges)

ref1 = (
    (np.sin(t[1:])**2 - np.sin(t[:-1])**2) / 2
) / dth

ref2 = (
    (np.cos(t[:-1])**3 - np.cos(t[1:])**3) / 3
) / dth

refphi = np.ones(N_PHI)

a1, c1, gl1, p1 = comparar(yth, eth, ref1)
a2, c2, gl2, p2 = comparar(yth, eth, ref2)
a3, c3, gl3, p3 = comparar(yph, eph, refphi)

print("\n=== VALIDACION ANGULAR CONJUNTA ===")
print("Cascadas:", N)
print("Muones seleccionados:", len(mu))
print()
print(f"Cenital cos(theta)sin(theta): "
      f"chi2={c1:.3f}, gl={gl1}, "
      f"chi2_red={c1/gl1:.3f}, p={p1:.6f}")
print(f"Cenital cos²(theta)sin(theta): "
      f"chi2={c2:.3f}, gl={gl2}, "
      f"chi2_red={c2/gl2:.3f}, p={p2:.6f}")
print(f"Azimut constante: "
      f"chi2={c3:.3f}, gl={gl3}, "
      f"chi2_red={c3/gl3:.3f}, p={p3:.6f}")

cth = (theta_edges[:-1] + theta_edges[1:]) / 2
cph = (phi_edges[:-1] + phi_edges[1:]) / 2

pd.DataFrame({
    "theta_min_deg": theta_edges[:-1],
    "theta_max_deg": theta_edges[1:],
    "theta_centro_deg": cth,
    "distribucion": yth,
    "incertidumbre": eth,
    "referencia_cos_sin": a1 * ref1,
    "referencia_cos2_sin": a2 * ref2
}).to_csv(SALIDA / "distribucion_cenital.csv", index=False)

pd.DataFrame({
    "phi_min_deg": phi_edges[:-1],
    "phi_max_deg": phi_edges[1:],
    "phi_centro_deg": cph,
    "distribucion": yph,
    "incertidumbre": eph,
    "referencia_constante": a3 * refphi
}).to_csv(SALIDA / "distribucion_azimutal.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.errorbar(cth, yth, yerr=eth, fmt="s",
            color="crimson", capsize=4,
            label="CORSIKA, 105.000 cascadas")
ax.plot(cth, a1 * ref1, "--", color="steelblue",
        label=r"$\cos\theta\sin\theta$")
ax.plot(cth, a2 * ref2, "-", color="black",
        label=r"$\cos^2\theta\sin\theta$")
ax.set(xlabel=r"Ángulo cenital $\theta$ (grados)",
       ylabel="Distribución ponderada (u. a./grado)",
       title=r"Muones al nivel del mar, $1\leq p\leq5$ GeV/c")
ax.grid(alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(SALIDA / "validacion_cenital.png", dpi=300)
fig.savefig(SALIDA / "validacion_cenital.pdf")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.errorbar(cph, yph, yerr=eph, fmt="s",
            color="crimson", capsize=4,
            label="CORSIKA, 105.000 cascadas")
ax.axhline(a3, color="black",
           label="Referencia constante")
ax.set(xlabel=r"Ángulo azimutal $\phi$ (grados)",
       ylabel="Distribución ponderada (u. a./grado)",
       title=r"Muones al nivel del mar, $1\leq p\leq5$ GeV/c")
ax.set_xlim(0, 360)
ax.grid(alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(SALIDA / "validacion_azimutal.png", dpi=300)
fig.savefig(SALIDA / "validacion_azimutal.pdf")
plt.close(fig)

print("\nGraficas guardadas en:", SALIDA)
