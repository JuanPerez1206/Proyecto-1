# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 10:49:13 2025

@author: JUAN PEREZ
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# entradas del usuario

lat = float(input("Latitud (ej. 6.27 para Medellín): "))
lon = float(input("Longitud (ej. -75.57 para Medellín): "))

año = int(input("Año (ej. 2025): "))
mes = int(input("Mes (ej. 3): "))
dia = int(input("Día (ej. 21): "))

tilt = float(input("Inclinación del panel en grados: "))
az_panel = float(input("Azimut del panel (180 = sur): "))

# rango horario del dia

zona = "America/Bogota"
fechas = pd.date_range(
    start=f"{año}-{mes:02d}-{dia:02d} 00:00:00",
    end=f"{año}-{mes:02d}-{dia:02d} 23:59:59",
    freq="1h", tz=zona
)

n_dia = fechas.dayofyear[0]


# posicion del sol

lat_rad = np.radians(lat)
tilt_rad = np.radians(tilt)
az_panel_rad = np.radians(az_panel)

decl = np.radians(23.45) * np.sin(np.radians(360*(284+n_dia)/365))

B = np.radians(360*(n_dia-81)/364)
EoT = 9.87*np.sin(2*B) - 7.53*np.cos(B) - 1.5*np.sin(B)

long_ref = -75
horas = fechas.hour + fechas.minute/60
AST = horas + (4*(long_ref - lon) + EoT)/60

H = np.radians(15*(AST-12))

altura = np.arcsin(np.sin(lat_rad)*np.sin(decl) + np.cos(lat_rad)*np.cos(decl)*np.cos(H))
zenit = np.pi/2 - altura

az_sol = np.arccos(
    (np.sin(decl)*np.cos(lat_rad) - np.cos(decl)*np.sin(lat_rad)*np.cos(H)) / np.cos(altura)
)
az_sol = np.where(H > 0, 2*np.pi - az_sol, az_sol)


# irradiancia extraterrestre

Gsc = 1367
E0 = 1 + 0.033*np.cos(2*np.pi*n_dia/365)
I0 = Gsc * E0

ghi = I0 * np.cos(zenit)
ghi = np.where(ghi > 0, ghi, 0)


# separar directa y difusa

dni = ghi * 0.8 / np.cos(zenit)
dni = np.where(ghi > 0, dni, 0)

dhi = ghi * 0.2


# irradiancia en el panel inclinado

cos_theta = (
    np.sin(decl)*np.sin(lat_rad)*np.cos(tilt_rad)
    - np.sin(decl)*np.cos(lat_rad)*np.sin(tilt_rad)*np.cos(az_panel_rad)
    + np.cos(decl)*np.cos(lat_rad)*np.cos(tilt_rad)*np.cos(H)
    + np.cos(decl)*np.sin(lat_rad)*np.sin(tilt_rad)*np.cos(az_panel_rad)*np.cos(H)
    + np.cos(decl)*np.sin(tilt_rad)*np.sin(az_panel_rad)*np.sin(H)
)

poa_dir = dni * np.maximum(cos_theta, 0)
poa_dif = dhi * (1+np.cos(tilt_rad))/2
poa_total = poa_dir + poa_dif


# factor climatico

condicion = input("\nCondición climática (soleado/nublado/lluvioso): ").strip().lower()

factores = {
    "soleado": 1.0,
    "nublado": 0.5,
    "lluvioso": 0.2
}

factor_clima = factores.get(condicion, 1.0)

# reducir irradiancia segun condicion climatica
poa_total = poa_total * factor_clima
poa_dir   = poa_dir   * factor_clima
poa_dif   = poa_dif   * factor_clima


# resultados

print("\nValores de irradiancia (05:30 - 20:00):")
tabla = pd.DataFrame({
    "Global": poa_total,
    "Directa": poa_dir,
    "Difusa": poa_dif
}, index=fechas)
print(tabla.between_time("05:30","20:00").head(10))

energia_Wh = np.trapz(poa_total, dx=1)
energia_kWh = energia_Wh/1000
print(f"\nEnergía solar del {año}-{mes:02d}-{dia:02d}:")
print(f"{energia_Wh:.2f} Wh/m²  ≈  {energia_kWh:.2f} kWh/m²")

# Gráfica

plt.figure(figsize=(11,5))
plt.plot(fechas, poa_total, label="Global (ajustada clima)")
plt.plot(fechas, poa_dir, label="Directa")
plt.plot(fechas, poa_dif, label="Difusa")
plt.xlabel("Hora")
plt.ylabel("Irradiancia (W/m²)")
plt.title(f"Irradiancia en el panel - {año}-{mes:02d}-{dia:02d}")
plt.legend()
plt.grid()
plt.xticks(rotation=45)
plt.show()

















