import argparse
import struct
import csv
import math
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Extraer muones de archivos binarios CORSIKA"
)
parser.add_argument("archivo", help="Ruta del archivo DAT")
parser.add_argument(
    "-o", "--output",
    help="Nombre del CSV de salida"
)
args = parser.parse_args()

archivo = Path(args.archivo)
salida = Path(args.output) if args.output else Path(
    f"muones_{archivo.name}.csv"
)

# Formato CORSIKA estándar sin thinning:
# 5733 floats por registro, 21 subbloques de 273 floats.
NWORDS = 5733
SUBBLOCK = 273
PARTICLE = 7

filas = []
shower_id = 0
n_eventos = 0
n_registros = 0

with archivo.open("rb") as f:
    while True:
        marcador = f.read(4)
        if not marcador:
            break

        if len(marcador) != 4:
            raise ValueError("Marcador Fortran incompleto")

        nbytes = struct.unpack("<i", marcador)[0]

        if nbytes != NWORDS * 4:
            raise ValueError(f"Tamaño de registro inesperado: {nbytes}")

        datos = f.read(nbytes)
        cierre = f.read(4)

        if len(datos) != nbytes or len(cierre) != 4:
            raise ValueError("Registro Fortran incompleto")

        if struct.unpack("<i", cierre)[0] != nbytes:
            raise ValueError("Marcadores Fortran inconsistentes")

        valores = struct.unpack("<5733f", datos)
        n_registros += 1

        for inicio in range(0, NWORDS, SUBBLOCK):
            bloque = valores[inicio:inicio + SUBBLOCK]
            identificador = bloque[0]

            if 217432 < identificador < 217434:
                shower_id += 1
                n_eventos += 1
                continue

            if (
                211284 < identificador < 211286
                or 52814 < identificador < 52816
                or 3396 < identificador < 3398
                or 3300 < identificador < 3302
            ):
                continue

            for j in range(0, SUBBLOCK, PARTICLE):
                p = bloque[j:j + PARTICLE]

                if len(p) != PARTICLE:
                    continue

                codigo = int(p[0])

                if codigo <= 0:
                    continue

                id_particula = codigo // 1000

                if id_particula not in (5, 6):
                    continue

                nivel = codigo % 10
                if nivel == 0:
                    nivel = 10

                px, py, pz = p[1], p[2], p[3]
                momento = math.sqrt(px**2 + py**2 + pz**2)

                if momento <= 0:
                    continue

                cos_theta = abs(pz) / momento
                theta = math.degrees(math.acos(
                    max(-1, min(1, cos_theta))
                ))
                phi = math.degrees(math.atan2(py, px))

                filas.append([
                    shower_id,
                    nivel,
                    id_particula,
                    px, py, pz,
                    momento,
                    cos_theta,
                    theta,
                    phi,
                    p[4], p[5], p[6]
                ])

columnas = [
    "shower_id", "nivel", "id_particula",
    "px", "py", "pz", "p",
    "cos_theta", "theta_deg", "phi_deg",
    "x_cm", "y_cm", "t_ns"
]

with salida.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(columnas)
    writer.writerows(filas)

print("=== EXTRACCION COMPLETADA ===")
print("Registros binarios:", n_registros)
print("Cascadas identificadas:", n_eventos)
print("Muones extraidos:", len(filas))

for nivel, altura in [(1, 450), (2, 30)]:
    seleccion = [r for r in filas if r[1] == nivel]
    positivos = sum(r[2] == 5 for r in seleccion)
    negativos = sum(r[2] == 6 for r in seleccion)

    print(f"\nNivel {nivel} ({altura} m)")
    print("  Muones positivos:", positivos)
    print("  Muones negativos:", negativos)
    print("  Total:", len(seleccion))

print("\nCSV generado:", salida)
