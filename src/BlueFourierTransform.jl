module BlueFourierTransform

using SpectraFromScratch
using BLUEs
using OffsetArrays

import SpectraFromScratch: FourierTransform

# Write your package code here.
# function blue_fourier_transform(
#     func, # takes input and get obs
#     y::Estimate, # observations
#     # x0::Estimate, # first guess
#     Ψ::FrequencySpectrum, # first guess frequency spectrum
#     σ2mean::Number) # uncertainty of mean value, not given by spectrum

function blue_fourier_transform(
    func, # takes input and get obs
    y, # observations
    Ψ, # first guess frequency spectrum
    σ2mean::Number) # uncertainty of mean value, not given by spectrum

    x0 = construct_first_guess(Ψ, σ2mean)
    
    # issue: what are the units?
    # maybe passing the first guess will solve it
    E = impulse_response(x0.v, func)
    println(E)
    
    # necessary?
    # first-guess prediction of the obs
    # y0 = E*x0

    # solve for Fourier Transform
    x1 = combine(x0, y, E)
    return x1
end

function construct_first_guess(Ψ, σ2mean)
    ffundamental = first(Ψ.freq) # fundamental frequency
    n = length(Ψ.psi)
    N = 2n # assume even number of points
    T = 1/ffundamental # repeat interval
    
    # A = amplitude of waves (?)
    # some care is taken for singular Nyquist frequency
    A0 = vcat(sqrt.(Ψ.psi[1:end-1] .* (N^2 / 2T)),
              sqrt.(Ψ.psi[end] .* (N^2 / T)))  

    sample_state_val = sqrt(first(Ψ.psi)*first(Ψ.freq))
    vx0 = zeros(eltype(sample_state_val),N)

    σx0 = vcat(N*σ2mean, #mean value (scaled to Fourier coefficient)
	A0, # real parts
	A0[1:end-1]) # imaginary parts w/o Nyquist

    x0 = Estimate(vx0, σx0)
end

function impulse_response(input0, func)
    output0 = vec(func(input0))
    E = zeros(eltype(first(output0)/first(input0)), length(output0), length(input0))
    input = deepcopy(input0)
    for i in 1:length(input0)
	delta = one(eltype(input))

	# delta = 1.0*unit(first(input))
	input[i] += delta
	E[:,i] = (vec(func(input)) -
		  vec(output0))./delta
	input[i] -= delta
    end
    return E
end

function SpectraFromScratch.FourierTransform(v::Vector, df::Number)
    N = length(v)
    n = SpectraFromScratch.fourier_modes(N) 

    # if iseven(N)
    #     n = -(N/2):((N/2)-1)
    # else
    #     n = -(N/2):((N/2))
    # end

    sample = first(v) + im*first(v)
    coeff = OffsetArray( fill( zero(eltype(sample)), N), n)

    N2 = Integer(N/2)
    
    # zero frequency coefficient is real
    coeff[0] = v[1]

    # Nyquist frequency coefficient is real
    coeff[-N2] = v[N2+1]

    # positive frequencies are complex
    coeff[1:N2-1] = v[2:N2] + im * v[(N2 + 2): end]

    # negative frequencies are a complex conjugate
    coeff[-1:-1:-N2+1] = v[2:N2] - im * v[(N2 +2):end]

    return FourierTransform(coeff, df)
end

end
