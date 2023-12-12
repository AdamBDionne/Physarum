#plotting
using PyCall
using PyPlot

function plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

    mpl = pyimport("matplotlib")
    rcParams = PyDict(matplotlib["rcParams"])
    rcParams["font.size"] = 18.0
    rcParams["font.family"] = "Avenir"
    rcParams["axes.linewidth"] = 2
    rcParams["figure.dpi"] = 300
    px = 1/plt.rcParams["figure.dpi"]

    mm = 1/25.4

    x_axis = [5*i/60 for i = 1:size(time_series_rms,2)]
    cmap = matplotlib.cm.get_cmap("plasma")

    #PLOT RMS MODEL
    fig = plt.figure(figsize = (2*2*45.75*mm, 2*2*30.875*mm))
    ax = plt.gca()
    ax1 = plt.subplot(221)
    ax3 = plt.subplot(223)
    ax2 = plt.subplot(222)
    ax4 = plt.subplot(224)

        for j = 1:Ne-3
            if j > Ne-11
                ax1.plot(x_axis, time_series_rms[j,:], color = cmap((Ne-j)/(Ne-1))[1:3], zorder=Ne-j)
            else
                ax1.plot(x_axis, time_series_rms[j,:], color = (cmap((Ne-j)/(Ne-1))[1:3]..., 0.25), zorder=Ne-j)
            end
        end

        for j = 1:Ne-3
            if j > Ne-11
                ax2.plot(x_axis, time_series_rms_model[j,:], color = cmap((Ne-j)/(Ne-1))[1:3], zorder=Ne-j)
            else
                ax2.plot(x_axis, time_series_rms_model[j,:], color = (cmap((Ne-j)/(Ne-1))[1:3]..., 0.25), zorder=Ne-j)
            end
        end


   
    j = argmax(avg_amp .* (1 .- (avg_amp .== maximum(avg_amp))))[1]
    ax1.plot(x_axis, time_series_rms[j,:], color = "C3", zorder=10000)
    j = argmax(avg_amp_model .* (1 .- (avg_amp_model .== maximum(avg_amp_model))))[1]
    ax2.plot(x_axis, time_series_rms_model[j,:], color = "C3", zorder=10000)
    
    j = argmax(avg_amp)[1]
    ax1.plot(x_axis, time_series_rms[j,:], color = "C2", zorder=10000)
    j = argmax(avg_amp_model)[1]
    ax2.plot(x_axis, time_series_rms_model[j,:], color = "C2",  zorder=10000)

  
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.yaxis.set_tick_params(length=6)
    ax2.xaxis.set_tick_params(width=2)
    ax2.yaxis.set_tick_params(width=2)
    ax2.xaxis.set_tick_params(length=6)
    ax2.yaxis.set_tick_params(length=6)

    ax1.xaxis.tick_top()
    ax2.xaxis.tick_top()
    ax2.yaxis.tick_right()
    ax4.yaxis.tick_right()



    for j = 1:Ne-3
        ax3.scatter(j,avg_amp[j],color = cmap((Ne-j)/(Ne-1))[1:3], s = 35)
        ax4.scatter(j,avg_amp_model[j],color = cmap((Ne-j)/(Ne-1))[1:3], s = 35)
    end


    j = argmax(avg_amp)[1]
    ax3.scatter(j,avg_amp[j], color = "C2", s = 50, marker = "D")
    j = argmax(avg_amp_model)[1]
    ax4.scatter(j,avg_amp_model[j], color = "C2", s = 50, marker = "D")
    j = argmax(avg_amp .* (1 .- (avg_amp .== maximum(avg_amp))))[1]
    ax3.scatter(j,avg_amp[j], color = "C3", s = 50, marker = "D")
    j = argmax(avg_amp_model .* (1 .- (avg_amp_model .== maximum(avg_amp_model))))[1]
    ax4.scatter(j,avg_amp_model[j], color = "C3", s = 50, marker = "D")
 

    ax3.xaxis.set_tick_params(width=2)
    ax3.yaxis.set_tick_params(width=2)
    ax3.xaxis.set_tick_params(length=6)
    ax3.yaxis.set_tick_params(length=6)
    ax4.xaxis.set_tick_params(width=2)
    ax4.yaxis.set_tick_params(width=2)
    ax4.xaxis.set_tick_params(length=6)
    ax4.yaxis.set_tick_params(length=6)


    #ax1.sharey(ax2)
    #ax4.sharey(ax3)

   #ax1.tick_params(labelleft=false)
    #ax3.tick_params(labelleft=false)
    #ax1.yaxis.set_tick_params(direction="in")
    #ax3.yaxis.set_tick_params(direction="in")

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    fig.supxlabel("Mode Number")
    fig.suptitle("Time (minutes)")
    ax1.set_ylabel("Projection \n RMS")
    ax3.set_ylabel("Projection \n Amplitude")

    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(top=0.8)
    
    ax1.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=2))
    ax2.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=2))
    ax3.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=2))
    ax4.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=2))

    plt.savefig(string("Modes/plots/",file_suffix,"_rms_model.pdf"), bbox_inches="tight")
end
