using Revise
using BlueFourierTransform
using Test
using BLUEs
using SpectraFromScratch
using OffsetArrays
using Unitful

@testset "BlueFourierTransform.jl" begin

    ## make a power law (sample frequency spectrum)
    ucases = (false, true)
    # case: units or no units
    for units in ucases
        
        # record length/ simulation length.

        if !units
            T = 1000 # repeat interval
            σ2 = 1.0  # estimate of total variance
            # contaminate observations with noise
            σn = 0.01
            σxbar = 1.0 # uncertainty of mean value, not given by spectrum
        else
            println("test with units")
            T = 1000u"yr" # repeat interval
            σ2 = 1.0u"K"^2  # estimate of total variance
            # contaminate observations with noise
            σn = 0.01u"K"
            σxbar = 1.0u"K" # uncertainty of mean value, not given by spectrum
        end

        N = 10 # number of points on regularly sampled grid
        β = 2.0  # power law 
        Δt = T/N # time interval of regular sampling
        Δf = Δt/T # Δf = bandwidth
        ffundamental = 1/T # fundamental frequency
        fnyquist = 1/2Δt # Nyquist frequency
        n = SpectraFromScratch.fourier_modes(N)
        nf = maximum(abs.(n))
        f_positive = (1:nf)/T #  goes between fundamental and Nyquist frequencies
        Ψ = spectral_power_law(f_positive, β, σ2) # power-law frequency spectrum

        ## make regular timeseries
        x̂true = SpectraFromScratch.FourierTransform(Ψ)
        xtrue = SpectraFromScratch.RegularTimeseries(x̂true)
        @test isapprox(sum(xtrue.x), zero(eltype(xtrue.x)),
            atol=1e-10*oneunit(eltype(xtrue.x)))
        
        ## will need to independently decide on the mean value
        # (i.e., no info in frequency spectrum)
        # here the mean is zero
    
        ## make irregular samples (using linear interpolation)
        # T_inclusive =  maximum(xtrue.time) - minimum(xtrue.time)
        M = 5; # number of observations
        t = rand(M)*T
    
        ## save samples as "perfect observations"
        # save the function that makes samples (linear interpolation)
        irregular_sample_fourier_transform(t) = [expand(t[i],x̂true) for i in eachindex(t)]
        ytrue = irregular_sample_fourier_transform(t)

        # also need the function for control variables (later move into master function)
        irregular_sample_control_variables(x) = [expand(t[i],FourierTransform(x,ffundamental)) for i in eachindex(t)]
        # y0 = irregular_sample_control_variables(u0.v) # interesting test but wrong scope to keep it
    
        n = σn * randn(M)
        y_contaminated = ytrue .+ n
    
        # save samples and expected noise together in an `Estimate`
        y = Estimate(y_contaminated, fill(σn, M))

        # lets try to split everything up
        x0 = BlueFourierTransform.construct_first_guess(Ψ, σxbar)

        E = BlueFourierTransform.impulse_response(x0.v,
            irregular_sample_control_variables)
    
        # even split up BLUEs.combine
        # x1 = combine(x0, y, E)

        nvec = y.v - E * x0.v
        sumP = y.P + E * x0.P * transpose(E)

        # needs AlgebraicArrays to be added in order to work
        v = x0.v + x0.P * transpose(E) *
            (  sumP \ nvec )
        P = x0.P - x0.P * transpose(E) *
            ( sumP \ (E * x0.P) )

        # Here, we solve everything at once
        # solve for BLUE of Fourier Transform
        u = FourierTransform(
            irregular_sample_control_variables, # takes input and get obs
            y, # observations
            Ψ, # first guess frequency spectrum
            σxbar) # uncertainty of mean value, not given by spectrum

        # how well does it recontruct the obs?
        ỹ = irregular_sample_control_variables(u.v)
        @test sqrt(sum(((y_contaminated-ỹ)/σn).^2)/M) < 1

        # how close is the Fourier Transform description of the data?
        # they actually look pretty different
        x̂̃ = FourierTransform(u.v, ffundamental)
        Δx_real = real.(x̂true.coeff) .- real.(x̂̃.coeff)
        Δx_imag = imag.(x̂true.coeff) .- imag.(x̂̃.coeff)

        # is it within the error bars? Would have to propagate the errors through `FourierTransform`
    
        # how faithful is it to the input spectrum?
        Ψ̃ = periodogram(x̂̃)    
    end
end
