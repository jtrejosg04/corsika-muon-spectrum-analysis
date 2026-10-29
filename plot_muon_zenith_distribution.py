from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
CARPETA = BASE / "resultados_angulares_conjunto"

df = pd.read_csv(CARPETA / "distribucion_cenital.csv")

fig, ax = plt.subplots(figsize=(8, 5.5))

# Resultados de CORSIKA con incertidumbres
ax.errorbar(
    df["theta_centro_deg"],
    df["distribucion"],
    yerr=df["incertidumbre"],
    fmt="s",
    color="crimson",
    capsize=4,
    label="CORSIKA, 105.000 cascadas"
)

# Unicamente la referencia cos²(theta) sin(theta)
ax.plot(
    df["theta_centro_deg"],
    df["referencia_cos2_sin"],
    "-",
    color="black",
    linewidth=1.8,
    label=r"$\cos^2\theta\sin\theta$"
)

ax.set_xlabel(r"Ángulo cenital $\theta$ (grados)")
ax.set_ylabel("Distribución ponderada (u. a./grado)")
ax.set_title(
    r"Muones al nivel del mar, $1\leq p\leq5$ GeV/c"
)

ax.grid(alpha=0.25)
ax.legend()

fig.tight_layout()

fig.savefig(
    CARPETA / "validacion_cenital_cos2.png",
    dpi=300
)
fig.savefig(
    CARPETA / "validacion_cenital_cos2.pdf"
)

plt.close(fig)

print("Gráfica guardada correctamente.")
