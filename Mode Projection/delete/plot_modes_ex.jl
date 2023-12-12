#plotting
using PyCall
using PyPlot
mpl = pyimport("matplotlib")

rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 18.0
rcParams["font.family"] = "Avenir"
rcParams["axes.linewidth"] = 2
rcParams["figure.dpi"] = 300
px = 1/plt.rcParams["figure.dpi"]

mm = 1/25.4

fig = plt.figure(figsize = (2*110*mm, 2*70*mm))
ax = plt.gca()

modes = copy(elastic_modes)
for j = 1:Ne
    modes[:,j] = elastic_modes[:,j] .- mean(elastic_modes[:,j])
    modes[:,j] = modes[:,j] ./ maximum(abs.(modes[:,j]))  
end
x_axis = [i for i = 1:Ne]

p1, = plt.plot(x_axis, modes[:,Ne-1], label="Mode N-1")
p2, = plt.plot(x_axis, modes[:,Ne-2], label="Mode N-2")
p3, = plt.plot(x_axis, modes[:,Ne-3], label="Mode N-3")


ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
mode="expand", borderaxespad=0, ncol=3,handles=[p1, p2, p3])
        # x_vals = [i for i = 1:Ne]
        # plt.plot(x_vals, odeSol[21][1:Ne])#, "#fc4a2b")
        # plt.plot(x_vals, odeSol[41][1:Ne])
        # plt.plot(x_vals, odeSol[61][1:Ne])
        # plt.plot(x_vals, odeSol[81][1:Ne])
plt.ylabel("Relative Amplitude (Au)")
plt.xlabel("Edge Index")
ax = plt.gca()
        # #ax.set_ylim([0, 1.05])

    ax.set_axisbelow(true)
    # alpha is used to control the transparency
    plt.grid(true, color="#93a1a1", alpha=0.3)
    
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(length=6)
    ax.yaxis.set_tick_params(length=6)
    
    
    plt.tight_layout()
    plt.savefig("Processed Data/path_graph/path_graph_modes.pdf")#, format="eps")