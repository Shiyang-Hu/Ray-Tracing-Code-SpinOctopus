import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm


x, y, z = np.loadtxt("data/j-profile-roz.dat", unpack=True)

x_min, x_max = -50, 50
y_min, y_max = -50, 50
radius = 2.0
gamma_value = 0.5
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
eht_cmap.set_bad(color='black')

plt.figure(figsize=(10, 10), facecolor='black')
ax = plt.gca()
ax.set_facecolor('black')


z_vmax = np.nanmax(z)

norm = PowerNorm(gamma=gamma_value, vmin=0, vmax=z_vmax)

plt.tricontourf(x, y, z, 
                levels=100, 
                cmap=eht_cmap, 
                norm=norm)

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

circle = plt.Circle((0, 0), radius, color='black', linestyle='-', linewidth=2, fill=True)
ax.add_artist(circle)

circle_edge = plt.Circle((0, 0), radius, color='white', linestyle='--', linewidth=2, fill=False)
ax.add_artist(circle_edge)

plt.gca().set_aspect('equal')
ax.set_xlabel(r'$r\cos\theta$', color='white', fontsize=20)
ax.set_ylabel(r'$r\sin\theta$', color='white', fontsize=20)

ax.tick_params(axis='both', which='major', colors='white', labelsize=12, width=1.5, length=6)
ax.tick_params(axis='both', which='minor', colors='white', labelsize=10, width=1, length=4)

ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white') 
ax.spines['left'].set_color('white')
ax.spines['right'].set_color('white')

#plt.text(1, 12.5, '(f)', fontsize=30, fontname='Times New Roman', color='white')

plt.gcf().set_facecolor('black')
plt.tight_layout()
plt.savefig('results/jet_emission.png', dpi=600, bbox_inches='tight', facecolor='black')
plt.show()