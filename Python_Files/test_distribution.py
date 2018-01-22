# -*- coding: utf-8 -*-
"""
test_distribution.py: module for testing kdes with probability density functions
"""

import numpy as np

from popsim import pdf
import matplotlib.pyplot as plt

xvals = np.linspace(0, 2*np.pi, 100)

theta = np.array([1,2,5,5.5])
Np = len(theta)
kappa = 30



thata_kde = pdf.estimate(theta, method="kde", varax=xvals, h=0.002, points=100)
thata_vM = pdf.estimate(theta, method="vonMises", varax=xvals, h=kappa, points=100)


plt.figure(figsize=(8,7))
plt.plot(thata_kde.varax,thata_kde.pdf)
plt.plot(thata_vM.varax,thata_vM.pdf)



for mu in theta:
    thata_vM = pdf.estimate(mu, method="vonMises", varax=xvals, h=kappa, points=100)
    plt.plot(thata_vM.varax,thata_vM.pdf/Np,'r--')
    

plt.show()
