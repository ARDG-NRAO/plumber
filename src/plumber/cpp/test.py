import ZernikePolyCalc_xt2py as zc
import numpy as np

c = np.linspace(0, 1, 67)
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)
#Z = zc.ZernikePolyCalc_xt2py(c, X, Y)
#print(Z)