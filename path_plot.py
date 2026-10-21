import numpy as np
import matplotlib.pyplot as plt

x, y, z = np.loadtxt("data/particle-path.dat", unpack=True)

fig1 = plt.figure(figsize=(10, 10))
ax = fig1.add_subplot(111, projection='3d')

R = 1.34117
spin = 0.94
u = np.linspace(0, 2 * np.pi, 180)
v = np.linspace(0, np.pi, 90)
xs = np.sqrt(R*R + spin*spin) * np.outer(np.cos(u), np.sin(v))
ys = np.sqrt(R*R + spin*spin) * np.outer(np.sin(u), np.sin(v))
zs = R * np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(xs, ys, zs, color='black', alpha=0.5, shade=True,
                linewidth=0, rstride=2, cstride=2)

ax.plot(x, y, z, color='blue', linewidth=1.5, alpha=0.8)
ax.scatter(x[0], y[0], z[0], color='green', s=50, marker='o', edgecolors='black')
ax.scatter(x[-1], y[-1], z[-1], color='red', s=50, marker='o', edgecolors='black')

ax.set_xlabel(r"$x^{\prime}$ [M]", fontsize=12, labelpad=10)
ax.set_ylabel(r"$y^{\prime}$ [M]", fontsize=12, labelpad=10)
ax.set_zlabel(r"$z^{\prime}$ [M]", fontsize=12, labelpad=10)

xmin, xmax = -20, 20
ymin, ymax = -20, 20
zmin, zmax = -20, 20
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
ax.set_box_aspect((xmax - xmin, ymax - ymin, zmax - zmin))

ax.view_init(elev=20, azim=-120)
ax.grid(True, alpha=0.3)
ax.text2D(0.02, 0.95, '(a)', transform=ax.transAxes,
          fontsize=40, fontname='Times New Roman', color='black')
plt.tight_layout()
plt.savefig('results/particle orbit.png', dpi=300, bbox_inches='tight')

n, hami_error = np.loadtxt("data/energy-error.dat", unpack=True)

fig2, ax2 = plt.subplots(figsize=(8, 5))
ax2.plot(n, hami_error, color='crimson', linewidth=1.5)

ax2.set_xlabel('n', fontsize=12)
ax2.set_ylabel('Hamiltonian error', fontsize=12)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('results/energy-error.png', dpi=300, bbox_inches='tight')

plt.show()