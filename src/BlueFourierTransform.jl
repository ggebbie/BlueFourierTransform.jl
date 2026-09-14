module BlueFourierTransform

using SpectraFromScratch
using BLUEs
    
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

    # necessary?
    # first-guess prediction of the obs
    # y0 = E*x0

    # solve for Fourier Transform
    x1 = combine(x0, y, E)
    return x1
end

function construct_first_guess(Ψ, σ2mean)
    ffundamental = first(Ψ.freq) # fundamental frequency
    N = 2length(Ψ.psi) # assume even number of points
    T = 1/ffundamental # repeat interval
    
    # A = amplitude of waves (?)
    A0 = sqrt.(Ψ.psi .* (N^2 / 2T))

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
		delta = 1.0*unit(first(input))
		input[i] += delta
		E[:,i] = (vec(func(input)) -
			vec(output0))./delta
		input[i] -= delta
	end
	return E
end

end
