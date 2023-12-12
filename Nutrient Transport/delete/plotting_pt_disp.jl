#plotting phase diagrams
using PyCall
using PyPlot
using DelimitedFiles
mpl = pyimport("matplotlib")
colors = pyimport("matplotlib.colors")

rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 18.0
rcParams["font.family"] = "Avenir"
rcParams["axes.linewidth"] = 2
rcParams["figure.dpi"] = 300
px = 1/plt.rcParams["figure.dpi"]

mm = 1/25.4

fig = plt.figure(figsize = (2*91.5*mm, 2*61.75*mm))
ax = plt.gca()



#get x and y axes
x_axis = [i for i = 1:20]
y_axis = [i for i = 1:20]

#color values
color_vals = t


p1 = plt.pcolormesh(x_axis, y_axis, color_vals .* Char_L .* 10^4, rasterized="true", cmap = "Spectral")

ax.xaxis.set_ticks([5.5, 10.5, 15.5])
ax.yaxis.set_ticks([5.5, 10.5, 15.5])
#ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
#ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
ax.set_xticklabels(["N-5", "N-10", "N-15"])
ax.set_yticklabels(["N-5", "N-10", "N-15"])
#ax.set_axisbelow(true)
# alpha is used to control the transparency
plt.grid(true, color="#93a1a1", alpha=0.3, zorder = 1000)

ax.xaxis.set_tick_params(width=2)
ax.yaxis.set_tick_params(width=2)
ax.xaxis.set_tick_params(length=6)
ax.yaxis.set_tick_params(length=6)

plt.xlabel("Mode Number")
plt.ylabel("Mode Number")


cbar = plt.colorbar(p1) #location ='bottom")
cbar.ax.set_ylabel("Average Displacement (mm)")
cbar.ax.tick_params(width = 2)
cbar.ax.tick_params(length = 6)
#plt.clim(minimum(t),0.9)
#cbar.ax.set_yticklabels(["0.7","0.8","0.9"])


#plt.colorbar(extend="max")
#plt.clim(0, 50);

plt.tight_layout()
plt.savefig("disp_pt_new.pdf")







