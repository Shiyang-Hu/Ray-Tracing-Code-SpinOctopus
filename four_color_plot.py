import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt("data/four-color-image.dat", unpack=True)
x_min, x_max = -15, 15
y_min, y_max = -15, 15
resolution = 200
color_mapping = {
    0: 'black',
    1: 'LawnGreen',
    2: 'yellow',
    3: 'DodgerBlue',
    4: 'OrangeRed',
    5: 'orange',
}
colors = [color_mapping[int(val)] for val in z]
num_x_bins = resolution
num_y_bins = resolution    
heatmap_data, xedges, yedges = np.histogram2d(x, y, bins=(num_x_bins, num_y_bins), 
                                               range=[[x.min(), x.max()], [y.min(), y.max()]], 
                                               weights=z)
color_indices = [list(color_mapping.keys()).index(int(val)) for val in z]
heatmap_colors = np.zeros_like(heatmap_data, dtype=int)
for i in range(len(color_indices)):
    x_idx = np.searchsorted(xedges, x[i]) - 1
    y_idx = np.searchsorted(yedges, y[i]) - 1
    heatmap_colors[x_idx, y_idx] = color_indices[i]
cmap = plt.cm.colors.ListedColormap(list(color_mapping.values()))
plt.figure(figsize=(10, 10))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.imshow(heatmap_colors.T, extent=extent, origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=len(color_mapping) - 1, interpolation='none')
plt.gca().set_aspect('equal')
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)
plt.gca().xaxis.set_tick_params(direction='in', top=True, bottom=True)
plt.gca().yaxis.set_tick_params(direction='in', left=True, right=True)
#plt.xlabel('x [M]', fontsize=18,fontname='Times New Roman')
#plt.ylabel('y [M]', fontsize=18,fontname='Times New Roman')
#plt.xticks(fontsize=14,fontname='Times New Roman')
#plt.yticks(fontsize=14,fontname='Times New Roman')
plt.text(-9, 7.5, '(c)', fontsize=40, fontname='Times New Roman', color='black')
plt.xticks([])
plt.yticks([])
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.savefig('results/four color image.png', dpi=600, bbox_inches='tight')
plt.show()





