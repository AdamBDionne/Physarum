
#file_suffix = "sample_spring_nine"
#plotting phase diagrams
using PyCall
using PyPlot
using DelimitedFiles
using Polynomials

function plot_disp(file_suffix)
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
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    #get data
    model_avg = readdlm(string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix,"_avg_disp.txt"))
    real_avg = readdlm(string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix,"_avg_disp.txt"))

    #get x and y axes
    x_axis_real = [6*i/60 for i = 1:size(real_avg,1)]
    x_axis_model = [6*i/60 for i = 1:size(model_avg,1)]

    model_avg = Char_L .* 10^4 .* model_avg[70:end]; real_avg = Char_L .* 10^4 .* real_avg[35:end]

    x_axis_real = x_axis_real[35:end]
    x_axis_model = x_axis_model[70:end]

    lg_model_avg = log.(model_avg); lg_real_avg = log.(real_avg);
    lg_x_axis_real = log.(x_axis_real); lg_x_axis_model = log.(x_axis_model)

    #find line of best fit for both 
        # fits = []
        # indices = []
        # for i = 150:20:size(model_avg,1)
        #     curr_fit = Polynomials.fit(lg_x_axis_model[1:i],lg_model_avg[1:i],1)
        #     push!(fits, curr_fit[1])
        #     push!(indices, i)
        # end
        # fit_until = argmax(fits)
        # fit_until = indices[fit_until]
    fit_until = size(model_avg,1)
    fit = Polynomials.fit(lg_x_axis_model[1:fit_until],lg_model_avg[1:fit_until],1)

    label_t = string("Slope = ",round(fit[1],digits=2))
    p1 = ax1.loglog(x_axis_model, model_avg,color="C4")
    p2 = ax1.loglog(x_axis_model[1:fit_until], exp.(fit.(lg_x_axis_model[1:fit_until])),label=label_t, linewidth=3,color="C0")
    ax1.minorticks_off()

    ax1.legend(loc=2)
    ax1.set_xlabel("Time (minutes)")
    ax1.set_ylabel("Displacement (mm)")
    ax1.set_title("Modeled")

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


    #find line of best fit for both 
        # fits = []
        # indices = []
        # for i = 150:20:size(real_avg,1)
        #     curr_fit = Polynomials.fit(lg_x_axis_real[1:i],lg_real_avg[1:i],1)
        #     push!(fits, curr_fit[1])
        #     push!(indices, i)
        # end
        # fit_until = argmax(fits)
        # fit_until = indices[fit_until]
    fit_until = size(real_avg,1)
    fit = Polynomials.fit(lg_x_axis_real[1:fit_until],lg_real_avg[1:fit_until],1)

    label_t = string("Slope = ",round(fit[1],digits=2))
    p1 = ax2.loglog(x_axis_real, real_avg, color="C4")
    p2 = ax2.loglog(x_axis_real[1:fit_until], exp.(fit.(lg_x_axis_real[1:fit_until])),label=label_t, linewidth=3, color="C0")
    ax2.minorticks_off()

    ax2.legend(loc=2)
    ax2.set_xlabel("Time (minutes)")
    ax2.set_ylabel("Displacement (mm)")
    ax2.set_title("Observed")

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
    plt.savefig(string(file_suffix,"_disp.pdf"))
end