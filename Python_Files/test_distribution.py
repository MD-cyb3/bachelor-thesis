# -*- coding: utf-8 -*-
"""
test_distribution.py: module for testing kdes with probability density functions
"""

import numpy as np

from popsim import pdf
import matplotlib.pyplot as plt

xvals = np.linspace(0, 2*np.pi, 100)

theta = np.array([1,2,5,5.5])



thata_kde = pdf.estimate(theta, method="kde", varax=xvals, points=100)
thata_vM = pdf.estimate(theta, method="vonMises", varax=xvals, h=8, points=100)


plt.figure(figsize=(8,7))
plt.plot(thata_kde.varax,thata_kde.pdf)
plt.plot(thata_vM.varax,thata_vM.pdf)
plt.show()

