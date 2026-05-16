# -*- coding: utf-8 -*-
import io

import pandas as pd


def analizar_tiempos_ejecucion(datos_csv):
    # 1. Cargar los datos usando BytesIO para evitar el error de unicode en Python 2
    df = pd.read_csv(io.BytesIO(datos_csv.strip()))

    # 2. Sustituir 'Asincrono' por 'Ivan'
    df["Version"] = df["Version"].replace("Asincrono", "Ivan")

    # 3. Agrupar por medias y pivotar (Compatible con Python 2.7)
    tabla_comparativa = pd.pivot_table(
        df,
        index=["Total_Procs", "N_Groups"],
        columns="Version",
        values="Tiempo_s",
        aggfunc="mean",
    ).reset_index()

    # 4. Redondear los tiempos a 3 decimales
    tabla_comparativa["Ivan"] = tabla_comparativa["Ivan"].round(3)
    tabla_comparativa["Jon"] = tabla_comparativa["Jon"].round(3)

    # 5. Calcular la diferencia (Quién gana en cada caso)
    tabla_comparativa["Mejor_Version"] = tabla_comparativa.apply(
        lambda row: "Ivan" if row["Ivan"] < row["Jon"] else "Jon", axis=1
    )

    return tabla_comparativa


# Tus datos de entrada (ya no hace falta poner la 'u' delante)
datos = """Version,Total_Procs,P_GroupSize,N_Groups,Iteracion,Tiempo_s
Asincrono,17,16,1,1,256.868
Asincrono,17,16,1,2,280.074
Asincrono,17,16,1,3,295.206
Jon,17,16,1,1,322.312
Jon,17,16,1,2,322.122
Jon,17,16,1,3,322.807
Asincrono,33,16,2,1,129.873
Asincrono,33,16,2,2,119.97
Asincrono,33,16,2,3,132.31
Jon,33,16,2,1,167.85
Jon,33,16,2,2,187.12
Jon,33,16,2,3,179.382
Asincrono,49,16,3,1,81.009
Asincrono,49,16,3,2,82.165
Asincrono,49,16,3,3,100.81
Jon,49,16,3,1,117.098
Jon,49,16,3,2,100.26
Jon,49,16,3,3,114.888
Asincrono,65,16,4,1,71.185
Asincrono,65,16,4,2,69.932
Asincrono,65,16,4,3,67.301
Jon,65,16,4,1,74.825
Jon,65,16,4,2,76.944
Jon,65,16,4,3,77.827
Asincrono,81,16,5,1,66.546
Asincrono,81,16,5,2,73.194
Asincrono,81,16,5,3,51.657
Jon,81,16,5,1,62.12
Jon,81,16,5,2,62.463
Jon,81,16,5,3,60.243
Asincrono,97,16,6,1,44.9
Asincrono,97,16,6,2,56.784
Asincrono,97,16,6,3,46.751
Jon,97,16,6,1,49.709
Jon,97,16,6,2,49.647
Jon,97,16,6,3,60.339
Asincrono,113,16,7,1,66.243
Asincrono,113,16,7,2,69.47
Asincrono,113,16,7,3,83.977
Jon,113,16,7,1,61.001
Jon,113,16,7,2,58.454
Jon,113,16,7,3,52.825"""

# Ejecutar y mostrar
resultado = analizar_tiempos_ejecucion(datos)

# Usamos un print clásico de la tabla por si el entorno no soporta to_markdown
print(resultado.to_string(index=False))
