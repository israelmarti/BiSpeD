#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bisped - Binary Spectral Disentangling (v1.7)
Interfaz pública del paquete.
"""

import multiprocessing as mp
import gc
import sys

# Importar todas las funciones públicas desde el núcleo
from bisped_core import (
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
    vexplore
)

# Configurar spawn
try:
    mp.set_start_method('spawn', force=True)
except RuntimeError:
    pass

if __name__ == "__main__":
    print("====================================\nBinary Spectral Disentangling (v1.7)\n====================================\n")
    print("Available functions list:\n")
    print("\tfind2c\n\thselect\n\tonecomp\n\tqfitg\n\trvbina\n\tsetrvs\n\tspbina\n\tsplot\n\tuniform\n\tvgrid\n\tvexplore\n\n")
