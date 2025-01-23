import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np


dimension = 20
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(1.5, 3.5, (3.5-1.5)/dimension)
Y = np.arange(1.8, 4.4, (4.4-1.8)/dimension)
X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
from icecream import ic
filename="out"
from os.path import exists
file_exists = exists(filename)
if file_exists:
    np_inData = np.genfromtxt(filename)
else:
    import os,subprocess
    ic(filename)
    subprocess.call("pwd",shell=True)
    os._exit(0)
Z = np_inData.reshape((dimension,dimension))



plt.xlabel('Cl-C distance',{'fontsize':16,'color':'black'})    # 設定 x 軸標籤
plt.ylabel('I-C distance',{'fontsize':16,'color':'black'})  # 設定 y 軸標籤
# Plot the surface.




# Make data
n = 100
xs = np.linspace(0, 1, n)
ys = np.sin(xs * 6 * np.pi)
zs = np.cos(xs * 6 * np.pi)

# Plot
ax.plot(xs, ys, zs)



surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
