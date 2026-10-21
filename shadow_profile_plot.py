import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.interpolate import splprep, splev
from scipy.optimize import least_squares
x, y = np.loadtxt("data/shadow.dat", unpack=True)
x_min_data = np.min(x)
x_max_data = np.max(x)
y_min_data = np.min(y)
y_max_data = np.max(y)
print(f"xmin: {x_min_data:.6f}")
print(f"xmax: {x_max_data:.6f}")
print(f"ymin: {y_min_data:.6f}")
print(f"ymax: {y_max_data:.6f}")
x_min, x_max = -10, 10
y_min, y_max = -10, 10  
plt.figure(figsize=(10, 10))
ax = plt.gca()
plt.scatter(x, y, s=6, alpha=0.7, color='black', edgecolors='none')
points = np.column_stack((x, y))
hull = ConvexHull(points)
hull_points = points[hull.vertices]
tck, u = splprep([hull_points[:, 0], hull_points[:, 1]], s=0.05, per=True)
u_new = np.linspace(0, 1, 1000)
x_contour, y_contour = splev(u_new, tck)
plt.plot(x_contour, y_contour, 'red', linewidth=3, alpha=0.9)
def circle_residuals(params, x, y):
    x0, y0, r = params
    return np.sqrt((x - x0)**2 + (y - y0)**2) - r
center_x = np.mean(x_contour)
center_y = np.mean(y_contour)
initial_radius = np.mean(np.sqrt((x_contour - center_x)**2 + (y_contour - center_y)**2))
initial_guess = [center_x, center_y, initial_radius]
result = least_squares(circle_residuals, initial_guess, args=(x_contour, y_contour))
x0_fit, y0_fit, r_fit = result.x
distances = np.sqrt((x_contour - x0_fit)**2 + (y_contour - y0_fit)**2)
std_radius = np.std(distances)
print(f"Center fitting: ({x0_fit:.6f}, {y0_fit:.6f})")
print(f"Radius fitting: {r_fit:.6f}")
print(f"Standard deviation: {std_radius:.6f}")
theta = np.linspace(0, 2*np.pi, 100)
circle_x = x0_fit + r_fit * np.cos(theta)
circle_y = y0_fit + r_fit * np.sin(theta)
#plt.plot(circle_x, circle_y, 'blue', linestyle='--', linewidth=2, alpha=0.8, label='Fitted circle')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal') 
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.xlabel('x [M]', fontsize=25)
plt.ylabel('y [M]', fontsize=25)
#plt.text(-18, 16, '(a)', fontsize=50, fontname='Times New Roman', color='black')
ax.tick_params(axis='both', which='major', labelsize=12)
plt.savefig('results/shadow contour.png', dpi=300, bbox_inches='tight')
plt.show()