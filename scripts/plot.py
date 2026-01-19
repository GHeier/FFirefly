import firefly as fly
import numpy as np
import matplotlib.pyplot as plt

fly.load_theme()

path = "GXMG"
npts = 50
field = fly.Field_C("data/base_chi.h5")

kpts = fly.plot.get_path(path, field.domain, npts)
x = np.arange(0, len(kpts))
y = np.array(field(kpts))

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.plot(x, y.real, label=f"T = {file_base[i]}")
fly.plot.format_path(ax, path, npts)
ax.legend()
plt.show()
plt.savefig("plots/chi.png")
