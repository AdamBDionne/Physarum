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

curr_ind = 2 

#source
source = "Phase Diagrams v.2/data/300_"
destination = "Phase Diagrams v.2/plots/300_"

source = ""
destination = ""

#get x and y axes
x_axis = vec(readdlm(string(source,"xaxis.txt")))
y_axis = vec(readdlm(string(source,"yaxis.txt")))

#color values
    #color_vals = readdlm(string(source,"mat_",curr_ind,".txt"))
    color_vals = readdlm("disp_pt.txt")

# for k = 1:5
#     for i = 2:size(color_vals,1)-1
#         for j = 2:size(color_vals,2)-1
#             if isnan(color_vals[i,j]) && sum(isnan.(color_vals[i-1:i+1,j-1:j+1])) < 4
#                 color_vals[i,j] = 0;
#                 A = color_vals[i-1:i+1,j-1:j+1]
#                 A = A[.!isnan.(A)]
#                 color_vals[i,j] = A[1]
#             end
#         end
#     end
# end

#x_crop = 25
y_crop = 16

#x_crop_ind = Int(floor(((100-x_crop)/(100)) * size(x_axis,1)))
y_crop_ind= Int(floor(((100-y_crop)/(100)) * size(y_axis,1)))
#x_axis = x_axis[1:x_crop_ind]
y_axis = y_axis[1:y_crop_ind]
#color_vals = color_vals[1:y_crop_ind,1:x_crop_ind]
color_vals = color_vals[1:y_crop_ind,:]
Z = color_vals

for i = 1:size(color_vals,1)
    for j = 1:size(color_vals,2)
        if color_vals[i,j] < 0
            color_vals[i,j] = 0 
        end
    end
end
p1 = plt.pcolormesh(x_axis, y_axis, color_vals .* Char_L .* 10^4, rasterized="true", cmap = "magma")#, norm=matplotlib.colors.LogNorm())#, #shading = "nearest")#, cmap = "viridis", interpolation = "nearest")


#ax.set_axisbelow(true)
# alpha is used to control the transparency
plt.grid(true, color="#93a1a1", alpha=0.3, zorder = 1000)

ax.xaxis.set_tick_params(width=2)
ax.yaxis.set_tick_params(width=2)
ax.xaxis.set_tick_params(length=6)
ax.yaxis.set_tick_params(length=6)

plt.xlabel("Contractile Strain")
plt.ylabel("Production Strain")

#formatter = matplotlib.ticker.LogFormatter(10, labelOnlyBase=false) 

t = color_vals[.!isnan.(color_vals)]

cbar = plt.colorbar(p1) #location ='bottom")
cbar.ax.set_ylabel("Average Displacement (mm)")
cbar.ax.tick_params(width = 2)
cbar.ax.tick_params(length = 6)
#plt.clim(minimum(t),0.9)
#cbar.ax.set_yticklabels(["0.7","0.8","0.9"])


#plt.colorbar(extend="max")
#plt.clim(0, 50);

plt.tight_layout()
plt.savefig(string(destination, "_", curr_ind,".eps"))







