import struct
import csv
import math
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("archivo")
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

filas = []

with open(args.archivo, "rb") as f:
    while True:
        marcador = f.read(4)
        if not marcador:
            break

        nbytes = struct.unpack("<i", marcador)[0]
        datos = f.read(nbytes)
        cierre = f.read(4)

        if nbytes != 5733 * 4:
            raise ValueError("Formato binario inesperado")

        if len(datos) != nbytes or len(cierre) != 4:
            raise ValueError("Registro incompleto")

        if struct.unpack("<i", cierre)[0] != nbytes:
            raise ValueError("Marcadores inconsistentes")

        valores = struct.unpack("<5733f", datos)

        for i in range(0, 5733, 273):
            bloque = valores[i:i+273]

            if 217432 < bloque[0] < 217434:
                shower_id = int(round(bloque[1]))
                tipo = int(round(bloque[2]))
                energia = bloque[3]
                theta = math.degrees(bloque[10])
                phi = math.degrees(bloque[11])

                filas.append([
                    shower_id, tipo, energia, theta, phi
                ])

with open(args.output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "shower_id", "primary_id",
        "energy_GeV", "theta_deg", "phi_deg"
    ])
    writer.writerows(filas)

print("=== CATALOGO DE PRIMARIOS ===")
print("Archivo:", args.output)
print("Primarios extraidos:", len(filas))
print("Primeros cinco registros:")

for fila in filas[:5]:
    print(fila)
