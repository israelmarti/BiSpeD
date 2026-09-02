#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bisped - Binary Spectral Disentangling (v1.7)
Interfaz pública del paquete.

Esta es la interfaz que los usuarios deben importar.
Todas las funciones principales están disponibles aquí.
"""

# Importar funciones públicas desde el núcleo
from .bisped_core import (
    find2c,
    hselect,
    onecomp,
    qfitg,
    rvbina,
    setrvs,
    spbina,
    splot,
    uniform,
    vgrid,
    vexplore,
)

# Versión del paquete
__version__ = "1.7"

# Banner de presentación (se muestra SOLO si se ejecuta el script directamente)
if __name__ == "__main__":
    print("====================================")
    print("Binary Spectral Disentangling (v1.7)")
    print("====================================\n")
    print("Available functions:\n")
    print("\tfind2c    - Cross-correlation for mass ratio and Teff")
    print("\thselect   - Extract FITS header keywords")
    print("\tonecomp   - Compare a spectrum with a template list")
    print("\tqfitg     - Gaussian fitting to q profiles")
    print("\trvbina    - Radial velocity iterative refinement")
    print("\tsetrvs    - Set radial velocities")
    print("\tspbina    - Compute primary and secondary spectra")
    print("\tsplot     - Plot a FITS spectrum")
    print("\tuniform   - Normalize spectra to a common continuum")
    print("\tvgrid     - Explore systemic RV grid")
    print("\tvexplore  - Explore vgrid results interactively")
