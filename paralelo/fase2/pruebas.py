# -*- coding: utf-8 -*-
from __future__ import print_function

import csv
import re
import subprocess
import sys

# Configuración de rutas (ajusta si los nombres de archivo son distintos)
EXEPAR = "../../../exepar"
CARD = "../card"

# Versiones a comparar
EJECUTABLES = {"Asincrono": "./asincrono/ejecutable", "Jon": "./jon/ejecutable"}

# Parámetros de las pruebas ajustados al "Punto Dulce" (P=16)
P_SIZES = [16]
# N es el número de grupos.
# Cálculos: 1+(1*16)=17, 1+(2*16)=33, ..., 1+(7*16)=113
N_GROUPS = [1, 2, 3, 4, 5, 6, 7]
REPETICIONES = 3  # Veces que se ejecuta cada prueba para hacer una media

OUTPUT_CSV = "comparativa_rendimiento_p16.csv"


def ejecutar_prueba(npr, exe_path, p_size):
    """Lanza el comando y extrae el tiempo de ejecución."""
    comando = ["bash", EXEPAR, str(npr), exe_path, CARD, str(p_size)]

    try:
        resultado_bytes = subprocess.check_output(comando, stderr=subprocess.STDOUT)
        resultado = resultado_bytes.decode("utf-8", errors="ignore")

        match = re.search(r"Time \(parallel\):\s+([\d\.]+)\s+s", resultado)
        if match:
            return float(match.group(1))
        else:
            print("Error: No se encontró el tiempo en la salida.")
            return None

    except subprocess.CalledProcessError as e:
        # ¡AQUÍ ESTÁ LA MAGIA! Si falla, imprimimos el mensaje real de OpenMPI
        mensaje_error = e.output.decode("utf-8", errors="ignore")
        print("\n[ERROR FATAL DE MPI]:\n" + mensaje_error)
        return None

    except Exception as e:
        print("Error inesperado: {}".format(e))
        return None


def main():
    print(
        "Iniciando pruebas de rendimiento (Punto Dulce P=16, Límite: 120 procesos)..."
    )

    with open(OUTPUT_CSV, mode="w") as file:
        writer = csv.writer(file)
        # Cabecera del CSV
        writer.writerow(
            [
                "Version",
                "Total_Procs",
                "P_GroupSize",
                "N_Groups",
                "Iteracion",
                "Tiempo_s",
            ]
        )

        for p in P_SIZES:
            for n in N_GROUPS:
                npr = 1 + (n * p)

                if npr > 120:
                    print("Saltando npr={} (excede el límite del clúster)".format(npr))
                    continue

                for nombre, path in EJECUTABLES.items():
                    # Usamos sys.stdout para que se imprima en la misma línea sin fallar en Python 2
                    mensaje = "Probando {} con {} procesos (P={}, N={})... ".format(
                        nombre, npr, p, n
                    )
                    sys.stdout.write(mensaje)
                    sys.stdout.flush()

                    for i in range(1, REPETICIONES + 1):
                        tiempo = ejecutar_prueba(npr, path, p)
                        if tiempo is not None:
                            writer.writerow([nombre, npr, p, n, i, tiempo])

                    print("Hecho.")

    print("\nPruebas finalizadas. Resultados guardados en: {}".format(OUTPUT_CSV))


if __name__ == "__main__":
    main()
