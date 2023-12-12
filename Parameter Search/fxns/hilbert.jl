function hilbert(x::AbstractArray{T,2}, n::Int64=0) where T <: Real
    x_ = copy(x)

    # work along columns
    size(x_,1)==1 ? x_ = permutedims(x_,[2 1]) : nothing

    if(n>0 && n<size(x_,1))
        x_ = x_[1:n,:]
    elseif(n>0 && n>size(x_,1))
        x_ = cat(1,x_,zeros(n-size(x_,1),size(x_,2)))
    else
        n = size(x_,1)
    end

    xf = fft(x_,1)
    h = zeros(Int64,n)
    if n>0 && n % 2 == 0
        #even, nonempty
        h[1:div(n,2)+1] .= 1
        h[2:div(n,2)] .= 2
    elseif n>0
        #odd, nonempty
        h[1] = 1
        h[2:div(n + 1,2)] .= 2
    end
    x_ = ifft(xf .* h[:,ones(Int64,size(xf,2))],1)

    # restore to original shape if necessary
    size(x,1)==1 ? x_ = permutedims(x_,[2 1]) : nothing

    return x_

end

function hilbert(x::AbstractArray{T}, n::Int64=0) where T <: Complex
    @warn("Using real part, ignoring complex part")
    Hilbert.hilbert(real(x))
end