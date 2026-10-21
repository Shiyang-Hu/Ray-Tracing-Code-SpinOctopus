import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.gridspec import GridSpec

x1, y1, z1 = np.loadtxt("data/j-profile-roz.dat", unpack=True)
x2, y2, z2 = np.loadtxt("data/j-profile-xoy.dat", unpack=True)
r,  f      = np.loadtxt("data/j-profile-r.dat",   unpack=True)

radius       = 1.34117
gamma_value  = 0.6
emission_max        = 1

x_min1, x_max1 = 0, 20
y_min1, y_max1 = -10, 10

x_min2, x_max2 = -20, 20
y_min2, y_max2 = -20, 20

r_min, r_max = 0, 30
f_min, f_max = 0, 1.5

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


def draw_profile(ax, x, y, z, xlim, ylim,
                 xlabel=r'$r\cos\theta$', ylabel=r'$r\sin\theta$'):
    ax.set_facecolor('black')
    norm = PowerNorm(gamma=gamma_value, vmin=0, vmax=emission_max)

    ax.tricontourf(x, y, z, levels=100, cmap=eht_cmap, norm=norm)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    ax.add_artist(plt.Circle((0, 0), radius, color='black',
                             linestyle='-', linewidth=2, fill=True))
    ax.add_artist(plt.Circle((0, 0), radius, color='white',
                             linestyle='--', linewidth=2, fill=False))

    ax.set_aspect('equal')
    ax.set_xlabel(xlabel, color='white', fontsize=18)
    ax.set_ylabel(ylabel, color='white', fontsize=18)

    ax.tick_params(axis='both', which='major', colors='white',
                   labelsize=11, width=1.5, length=6)
    ax.tick_params(axis='both', which='minor', colors='white',
                   labelsize=9, width=1, length=4)
    for sp in ax.spines.values():
        sp.set_color('white')

fig = plt.figure(figsize=(15, 11), facecolor='black')

gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])

draw_profile(ax1, x1, y1, z1,
             (x_min1, x_max1), (y_min1, y_max1),
             xlabel=r'$r\cos\theta$', ylabel=r'$r\sin\theta$')

draw_profile(ax2, x2, y2, z2,
             (x_min2, x_max2), (y_min2, y_max2),
             xlabel=r'$x$ [M]', ylabel=r'$y$ [M]')

ax3.set_facecolor('black')

ax3.plot(r, f, color='#FF6600', linewidth=2.0, zorder=2)
ax3.scatter(r, f, color='#FFCC00', s=12, zorder=3, edgecolors='none')

ax3.axvline(radius, color='white', linestyle='--',
            linewidth=1.2, alpha=0.7)

ax3.set_xlim(r_min, r_max)
ax3.set_ylim(f_min, f_max)
ax3.set_xlabel(r'$r$ [M]', color='white', fontsize=18)
ax3.set_ylabel(r'$f$',     color='white', fontsize=18)

ax3.tick_params(axis='both', which='major', colors='white',
                labelsize=11, width=1.5, length=6)
ax3.tick_params(axis='both', which='minor', colors='white',
                labelsize=9, width=1, length=4)
for sp in ax3.spines.values():
    sp.set_color('white')

ax3.grid(True, color='gray', alpha=0.3, linestyle=':')

plt.tight_layout(pad=0.4, h_pad=0.3, w_pad=0.3)

plt.savefig('results/disk_emission.png',
            dpi=600,
            bbox_inches='tight',
            pad_inches=0.05,
            facecolor='black')
plt.show()