import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

x, y, z = np.loadtxt("data/image.dat", unpack=True)
x_min, x_max = -15, 15
y_min, y_max = -15, 15
resolution = 200
gamma_value = 1
v_max=np.nanmax(z)*0.8
#v_max=9.0
x_edges = np.linspace(x_min, x_max, resolution)
y_edges = np.linspace(y_min, y_max, resolution)
X, Y = np.meshgrid(x_edges, y_edges)
Z_grid = np.full_like(X, np.nan, dtype=float)

z_vmax = np.nanmax(z)

for i, (x_val, y_val, z_val) in enumerate(zip(x, y, z)):
    x_index = np.argmin(np.abs(x_edges - x_val))
    y_index = np.argmin(np.abs(y_edges - y_val))
    Z_grid[y_index, x_index] = z_val

eht_colors = [
    (0.00, '#000000'),
    (0.15, '#330000'),  
    (0.30, '#660000'),  
    (0.45, '#990000'),  
    (0.60, '#CC0000'),  
    (0.75, '#FF3300'),  
    (0.85, '#FF6600'),  
    (0.92, '#FF9900'),  
    (0.97, '#FFCC00'),  
    (1.00, '#FFFFFF')   
]

eht_cmap = LinearSegmentedColormap.from_list('eht_black_hole', eht_colors)

norm = PowerNorm(gamma=gamma_value, vmin=0, vmax=v_max)

plt.figure(figsize=(10, 10))
Z_masked = np.ma.masked_invalid(Z_grid)
eht_cmap.set_bad(color='black')

plt.imshow(Z_masked, 
           cmap=eht_cmap, 
           norm=norm,  
           origin='lower', 
           extent=[x_min, x_max, y_min, y_max], 
           interpolation='none', 
           aspect='auto')

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal')
plt.axis('off')
#plt.gca().set_facecolor('black')
#plt.gcf().set_facecolor('black')
plt.tight_layout()
plt.text(-14, 12.5, '(c)', fontsize=40, fontname='Times New Roman', color='white')
print(f"Max Intensity: {z_vmax}")
print(f"Gamma value: {gamma_value}")
plt.savefig('results/image.png', dpi=600, bbox_inches='tight')
plt.show()