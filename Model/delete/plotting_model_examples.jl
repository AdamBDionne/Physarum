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

fig = plt.figure(figsize = (2*150*mm, 2*45*mm))
ax = plt.gca()

#fig = plt.figure(figsize=(2*750*px, 2*500*px), dpi=300)

        # x_vals = [i for i = 1:Ne]
        # plt.plot(x_vals, odeSol[21][1:Ne])#, "#fc4a2b")
        # plt.plot(x_vals, odeSol[41][1:Ne])
        # plt.plot(x_vals, odeSol[61][1:Ne])
        # plt.plot(x_vals, odeSol[81][1:Ne])
        # plt.ylabel("Volume (Au)")
        # ax = plt.gca()
        # #ax.set_ylim([0, 1.05])
        # plt.xlabel("Edge Index")
        # #highlight_num = argmax(local_order_signal[100:end-100])+100
        # #plt.plot(time_s[highlight_num], local_order_signal[highlight_num], marker="*", color="#0E5D11", markersize = 15)

        # ax.set_axisbelow(true)
        # # alpha is used to control the transparency
        # plt.grid(true, color="#93a1a1", alpha=0.3)

        # ax.xaxis.set_tick_params(width=2)
        # ax.yaxis.set_tick_params(width=2)
        # ax.xaxis.set_tick_params(length=6)
        # ax.yaxis.set_tick_params(length=6)


        # plt.tight_layout()
        # plt.savefig(string(destination, "_spatial_cord.pdf"))#, format="eps")


                # x_vals = [(i/fps)*15 for i in 1:fps*(t_f-t_i)]

                # signals = zeros(size(odeSol,1),2Ne)
                # for i = 1:size(odeSol,1)
                #     signals[i,:] = odeSol[i]
                # end

                # cmap = matplotlib.cm.get_cmap("Spectral")

                # for i = 1:Ne
                #     plt.plot(x_vals, signals[:,i], color = cmap(i/Ne))
                # end
                # plt.ylabel("Volume (Au)")
                # ax = plt.gca()
                # #ax.set_ylim([0, 1.05])
                # plt.xlabel("Time (s)")
                # #highlight_num = argmax(local_order_signal[100:end-100])+100
                # #plt.plot(time_s[highlight_num], local_order_signal[highlight_num], marker="*", color="#0E5D11", markersize = 15)

                # ax.set_axisbelow(true)
                # # alpha is used to control the transparency
                # plt.grid(true, color="#93a1a1", alpha=0.3)

                # ax.xaxis.set_tick_params(width=2)
                # ax.yaxis.set_tick_params(width=2)
                # ax.xaxis.set_tick_params(length=6)
                # ax.yaxis.set_tick_params(length=6)


                # plt.tight_layout()
                # plt.savefig(string(destination, "_vol_time.pdf"))#, format="eps")



    # x_vals = [(i/fps)*15 for i in 1:fps*(t_f-t_i)]

    # signals = zeros(size(odeSol,1),2Ne)
    # for i = 1:size(odeSol,1)
    #     signals[i,:] = odeSol[i]
    # end

    # cmap = matplotlib.cm.get_cmap("Spectral")

    # for i = 1:Ne
    #     plt.plot(x_vals, signals[:,Ne+i], color = cmap(i/Ne))
    # end
    # plt.ylabel("Concentration (Au)")
    # ax = plt.gca()
    # #ax.set_ylim([0, 1.05])
    # plt.xlabel("Time (s)")
    # #highlight_num = argmax(local_order_signal[100:end-100])+100
    # #plt.plot(time_s[highlight_num], local_order_signal[highlight_num], marker="*", color="#0E5D11", markersize = 15)

    # ax.set_axisbelow(true)
    # # alpha is used to control the transparency
    # plt.grid(true, color="#93a1a1", alpha=0.3)

    # ax.xaxis.set_tick_params(width=2)
    # ax.yaxis.set_tick_params(width=2)
    # ax.xaxis.set_tick_params(length=6)
    # ax.yaxis.set_tick_params(length=6)


    # plt.tight_layout()
    # plt.savefig(string(destination, "_conc_time.pdf"))#, format="eps")



    # x_vals = [(i/fps)*15 for i in 1:fps*(t_f-t_i)]

    # signals = zeros(size(odeSol,1),Ne)
    # for i = 1:size(odeSol,1)
    #     signals[i,:] = 2 * flows[i] ./ (pi .* avg_radii.^2)
    # end
    
    # cmap = matplotlib.cm.get_cmap("Spectral")
    
    # for i = 1:Ne
    #     plt.plot(x_vals, signals[:,i], color = cmap(i/Ne))
    # end
    # plt.ylabel("Flow Velocity (Au)")
    # ax = plt.gca()
    # #ax.set_ylim([0, 1.05])
    # plt.xlabel("Time (s)")
    # #highlight_num = argmax(local_order_signal[100:end-100])+100
    # #plt.plot(time_s[highlight_num], local_order_signal[highlight_num], marker="*", color="#0E5D11", markersize = 15)
    
    # ax.set_axisbelow(true)
    # # alpha is used to control the transparency
    # plt.grid(true, color="#93a1a1", alpha=0.3)
    
    # ax.xaxis.set_tick_params(width=2)
    # ax.yaxis.set_tick_params(width=2)
    # ax.xaxis.set_tick_params(length=6)
    # ax.yaxis.set_tick_params(length=6)
    
    
    # plt.tight_layout()
    # plt.savefig(string(destination, "_flow_time.pdf"))#, format="eps")
    
    
    numFrames = 900
    edge_ind = 52

    x_vals = [(i/fps)*15 for i in 1:numFrames]

    #get volume, conc, flow
    vol = zeros(numFrames,1)
    conc = zeros(numFrames,1)
    flow = zeros(numFrames,1)
    for i = 1:numFrames
        vol[i] = odeSol[i][edge_ind]
        conc[i] = odeSol[i][Ne+edge_ind]
        flow[i] = 2 * flows[i][edge_ind] / (pi * avg_radii[edge_ind].^2)
    end

    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin2.spines.right.set_position(("axes", 1.2))

    p1, = ax.plot(x_vals,vol,label="Volume")
    p2, = twin1.plot(x_vals,conc,"C1",label="Concentration")
    p3, = twin2.plot(x_vals,flow,"C2",label="Flow Velocity")
    
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Volume (Au)")
    twin1.set_ylabel("Concentration (Au)")
    twin2.set_ylabel("Flow Velocity (Au)")

    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    twin2.yaxis.set_ticks([4, 1, -2])

    # tkw = dict(size=4, width=1.5)
    # ax.tick_params(axis='y', colors=p1.get_color(), tkw)
    # twin1.tick_params(axis='y', colors=p2.get_color(), tkw)
    # twin2.tick_params(axis='y', colors=p3.get_color(), tkw)
    # ax.tick_params(axis='x', tkw)

    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
    mode="expand", borderaxespad=0, ncol=3,handles=[p1, p2, p3])

    #ax = plt.gca()
   

    ax.set_axisbelow(true)
    # alpha is used to control the transparency
    plt.grid(true, color="#93a1a1", alpha=0.3)
    
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(length=6)
    ax.yaxis.set_tick_params(length=6)
    
    twin1.xaxis.set_tick_params(width=2)
    twin1.yaxis.set_tick_params(width=2)
    twin1.xaxis.set_tick_params(length=6)
    twin1.yaxis.set_tick_params(length=6)

    twin2.xaxis.set_tick_params(width=2)
    twin2.yaxis.set_tick_params(width=2)
    twin2.xaxis.set_tick_params(length=6)
    twin2.yaxis.set_tick_params(length=6)
    
    plt.tight_layout()
    plt.savefig(string(destination, "_example_edge.pdf"))#, format="eps")
    
    
    