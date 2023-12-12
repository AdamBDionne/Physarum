using PyPlot
using PyCall

mpimg = pyimport("matplotlib.image")
# Plot results 

PyPlot.close()
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 24

fig, ax = PyPlot.subplots(figsize=(14,9))

num1 = Ne-4
num2 = Ne-5

#largest mode

fig.suptitle(string("Sample One, with Order Parameter ",round(mean(order_time_series),digits=3),"\n and Phase Difference between Two Dominant Harmonics ", round(phase_dif,digits=3)))

sub1 = plt.subplot(1,4,1)
img_mode= mpimg.imread(string("Processed Data/",file_suffix,"/",file_suffix,"mode_",num1,".png"))
sub1.imshow(img_mode)
sub1.axis("off")

#second largest mode 
sub2 = plt.subplot(1,4,2)
img_mode = mpimg.imread(string("Processed Data/",file_suffix,"/",file_suffix,"mode_",num2,".png"))
sub2.imshow(img_mode)
sub2.axis("off")

time = [5*i/60 for i = 1:numSamples]
# time series RMS
plt.subplot(1,4,3)
for i = 1:size(modes)[1]-1
    PyPlot.plot(time, abs.(time_series_rms[i,:]))
end
PyPlot.xlabel("Time (minutes)")
PyPlot.ylabel("RMS Mode Projection")

# organization 
plt.subplot(1,4,4)
for i = 1:size(modes)[1]-1
    PyPlot.scatter(sums_ordered[i,1], (sums_ordered[i,2]))
end
PyPlot.xlabel("Mode Number")
PyPlot.ylabel("Projection Size")

# title with phase difference & order parameter




plt.tight_layout()
PyPlot.savefig(string("mode_projections", file_suffix,".png"))