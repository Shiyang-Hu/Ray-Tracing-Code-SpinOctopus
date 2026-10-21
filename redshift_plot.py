import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt("data/redshift-1.dat", unpack=True)
x_min, x_max = -15, 15
y_min, y_max = -15, 15
resolution = 1000
x_edges = np.linspace(x_min, x_max, resolution)
y_edges = np.linspace(y_min, y_max, resolution)
X, Y = np.meshgrid(x_edges, y_edges)
Z_grid = np.full_like(X, np.nan, dtype=float)
for i, (x_val, y_val, z_val) in enumerate(zip(x, y, z)):
    x_index = np.argmin(np.abs(x_edges - x_val))
    y_index = np.argmin(np.abs(y_edges - y_val))
    Z_grid[y_index, x_index] = z_val
plt.figure(figsize=(10, 10))
cmap = plt.cm.coolwarm.reversed()
Z_masked = np.ma.masked_invalid(Z_grid)
cmap.set_bad(color='black')
plt.imshow(Z_masked, cmap=cmap, origin='lower', extent=[x_min, x_max, y_min, y_max], interpolation='none', vmin=0, vmax=2, aspect='auto')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal')
plt.axis('off')
plt.text(-14, 12.5, '(c)', fontsize=40, fontname='Times New Roman', color='black')
#plt.tight_layout()
plt.savefig('results/redshift.png', dpi=600, bbox_inches='tight')
plt.show()