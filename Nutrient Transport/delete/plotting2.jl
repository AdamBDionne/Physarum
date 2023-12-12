
file_suffix = "sample_spring_nine"

#plotting phase diagrams
using PyCall
using PyPlot
using DelimitedFiles
using Polynomials
mpl = pyimport("matplotlib")
colors = pyimport("matplotlib.colors")

rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 18.0
rcParams["font.family"] = "Avenir"
rcParams["axes.linewidth"] = 2
rcParams["figure.dpi"] = 300
px = 1/plt.rcParams["figure.dpi"]

mm = 1/25.4

fig = plt.figure(figsize = (4*91.5*mm, 2*61.75*mm))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

t_s = avg_disp;
avg_displacement = disp0;

#get x and y axes
x_axis_real = [10*i/60 for i = 1:size(avg_displacement,1)]
t_s = t_s .* Char_L .* 10^4

ex1 = disp0 .* Char_L .* 10^4
ex2 = disp1 .* Char_L .* 10^4
ex3 = disp2 .* Char_L .* 10^4
ex4 = disp3 .* Char_L .* 10^4

x_axis_model = [i for i = -pi:0.05:pi]
ax1.scatter(x_axis_model, t_s, color="C4")

ax1.scatter(x_axis_model[64], t_s[64], label="0 Phase Shift", color = "C0", s = 130, marker = "D")
ax1.scatter(x_axis_model[74], t_s[74], label="π/6 Phase Shift", color = "C1", s = 130, marker = "D")
ax1.scatter(x_axis_model[85], t_s[85], label="π/3 Phase Shift", color = "C2", s = 130, marker = "D")
ax1.scatter(x_axis_model[95], t_s[95], label="π/2 Phase Shift", color = "C3", s = 130, marker = "D")



ax1.set_xlabel("Phase (radians)")
ax1.set_ylabel("Displacement (mm)")
#ax1.set_title("Two Most Dominant Modes", y = 1.22)
ax1.legend(bbox_to_anchor=(0,1.1,1,0.2), loc="center", borderaxespad=0, ncol=2)

ax1.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3))
ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3))
ax1.grid(true, color="#93a1a1", alpha=0.3, zorder = 1000)

ax1.xaxis.set_tick_params(width=2)
ax1.yaxis.set_tick_params(width=2)
ax1.xaxis.set_tick_params(length=6)
ax1.yaxis.set_tick_params(length=6)

ax1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax1.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax1.set_xticks([-pi, 0, pi])
labels = ["-π", "0", "π"]
ax1.set_xticklabels(labels)

x_axis_real = x_axis_real[35:end]
ex1 = ex1[35:end]
ex2 = ex2[35:end]
ex3 = ex3[35:end]
ex4 = ex4[35:end]

lg_x_axis_real = log.(x_axis_real)
lg_ex1 = log.(ex1)
lg_ex2 = log.(ex2)
lg_ex3 = log.(ex3)
lg_ex4 = log.(ex4)

#find line of best fit for all
fit1 = Polynomials.fit(lg_x_axis_real,lg_ex1,1)
fit2 = Polynomials.fit(lg_x_axis_real,lg_ex2,1)
fit3 = Polynomials.fit(lg_x_axis_real,lg_ex3,1)
fit4 = Polynomials.fit(lg_x_axis_real,lg_ex4,1)

label1 = string("Slope = ",round(fit1[1],digits=2))
label2 = string("Slope = ",round(fit2[1],digits=2))
label3 = string("Slope = ",round(fit3[1],digits=2))
label4 = string("Slope = ",round(fit4[1],digits=2))

ax2.loglog(x_axis_real, ex1, color="C0",label=label1)
ax2.loglog(x_axis_real, exp.(fit1.(lg_x_axis_real)), linewidth=2, color="C6",alpha=0.5)

ax2.loglog(x_axis_real, ex2, color="C1",label=label2)
ax2.loglog(x_axis_real, exp.(fit2.(lg_x_axis_real)), linewidth=2, color="C6",alpha=0.5)

ax2.loglog(x_axis_real, ex3, color="C2",label=label3)
ax2.loglog(x_axis_real, exp.(fit3.(lg_x_axis_real)), linewidth=2, color="C6", alpha=0.5)

ax2.loglog(x_axis_real, ex4, color="C3",label=label4)
ax2.loglog(x_axis_real, exp.(fit4.(lg_x_axis_real)), linewidth=2, color="C6", alpha=0.5)



ax2.minorticks_off()

ax2.legend(bbox_to_anchor=(0,1.1,1,0.2), loc="center", borderaxespad=0, ncol=2)
ax2.set_xlabel("Time (minutes)")
ax2.set_ylabel("Displacement (mm)")
#ax2.set_title("Phase Shift = π/2", y = 1.22)

ax2.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3))
ax2.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3))
ax2.grid(true, color="#93a1a1", alpha=0.3, zorder = 1000)

ax2.xaxis.set_tick_params(width=2)
ax2.yaxis.set_tick_params(width=2)
ax2.xaxis.set_tick_params(length=6)
ax2.yaxis.set_tick_params(length=6)

ax2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax2.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

plt.tight_layout()
plt.savefig(string("stuff_disp.pdf"))