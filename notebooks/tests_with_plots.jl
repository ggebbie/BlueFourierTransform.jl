### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 8b2fe3cb-4557-4634-82b8-1e14e403bfcb
import Pkg; Pkg.activate(".")

# ╔═╡ 8cbfce62-ea81-42ed-95bd-3bfc052c50ac
Pkg.develop(path="../../SpectraFromScratch.jl")

# ╔═╡ 09001ed8-00a9-4d65-bf2f-70a87fac3e4f
begin
	using Revise
	using BlueFourierTransform
	using Test
	using BLUEs
	using SpectraFromScratch
	using OffsetArrays
	using Plots
	using Interpolations
end

# ╔═╡ 4b0a896e-b28a-11f1-b73f-77977d857f62
# testing BlueFourierTransform

# ╔═╡ 1e8fe281-32a0-4923-8b1a-b35d6bdc28a7
plotly()

# ╔═╡ 5ca1234f-374e-4838-834d-96a792125b9a
begin
    ## make a power law (sample frequency spectrum)

    # record length/ simulation length.
    T = 1000 # repeat interval
    N = 20 # number of points on regularly sampled grid
    β = 2.0  # power law 
    Δt = T/N # time interval of regular sampling
    fwidth = N/T # Δf = bandwidth, min to max frequencies
    ffundamental = 1/T # fundamental frequency
    fnyquist = 1/2Δt # Nyquist frequency
    σ2 = 1.0  # estimate of total variance
    n = SpectraFromScratch.fourier_modes(N)
    nf = maximum(abs.(n))
    f_positive = (1:nf)/T #  goes between fundamental and Nyquist frequencies
    Ψ = spectral_power_law(f_positive, β, σ2) # power-law frequency spectrum
end

# ╔═╡ 16f76114-2812-4ee0-b6bf-6b9832d08757
# plot spectrum
plot(Ψ.freq,
	Ψ.psi,
	#ylims=(1e-5,1e3),
	yscale=:log10,
 	xscale=:log10,
	color=:black,
	lw=4,
	title="frequency spectra",
	label="power law"
	)

# ╔═╡ f7a1395a-b65a-4d84-a062-8851b82af3ad
## make regular timeseries (non-deterministic)
x̂true = SpectraFromScratch.FourierTransform(Ψ)    

# ╔═╡ 81ea80de-4e29-46e5-a5b2-3ab12faee0b5
xtrue = SpectraFromScratch.RegularTimeseries(x̂true)

# ╔═╡ 4d278660-dd7a-4749-9765-67ad369eefa2
plot(0.0 .+ xtrue.time,
	xtrue.x,
	title = "true timeseries",
	legend=false)

# ╔═╡ 7ce0b464-a849-4211-8a14-314fb4294eb3
## make irregular samples (using linear interpolation)
begin
    Tinclusive =  maximum(xtrue.time) - minimum(xtrue.time)
    M = 1000; # number of observations
    t = rand(M)*Tinclusive
end

# ╔═╡ 6e0e263e-4070-4d00-9769-769c03bc2ef2
## save samples as "perfect observations"
# save the function that makes samples (linear interpolation)
begin
    irregular_sample_fourier_transform(t) = [expand(t[i],x̂true) for i in eachindex(t)]
    ytrue = irregular_sample_fourier_transform(t)
end

# ╔═╡ dff337bb-9954-4663-8b3f-d34cc16a44b6
# also need the function for control variables (later move into master function)
irregular_sample_control_variables(u) = 
    [expand(t[i],FourierTransform(u,ffundamental)) for i in eachindex(t)]

# ╔═╡ 32724079-1a28-4aa0-8dba-2851471e38ed
# lo-res way to get variable at any time
function irregular_sample_fourier_transform_interp(t, beta::FourierTransform)
    x = RegularTimeseries(beta)
    itp = linear_interpolation(
            OffsetArrays.no_offset_view(x.time),
            OffsetArrays.no_offset_view(x.x))
    out = [itp(i) for i in t]
end

# ╔═╡ bcd66b85-4c58-41d5-8958-b79142080f1d
ytrue_interp = irregular_sample_fourier_transform_interp(t, x̂true)

# ╔═╡ 8b66a95a-1e69-4ade-9eec-da0a1d7aef2c
begin
    # contaminate observations with noise
    σn = 0.01
    noise = σn * randn(M)
    y = ytrue .+ noise
    y_interp = ytrue_interp .+ noise
end

# ╔═╡ c746ff66-9c81-48d7-86a8-aee1a2c05dc8
begin
	plot(0.0 .+ xtrue.time,
	xtrue.x,
		 label = "linear interp. truth",
	title = "timeseries and obs",
		 lw = 5,
	legend=true)

    Plots.scatter!(t,
					y,
					yerror=σn,
					label="obs")	
	Plots.scatter!(t,
					y_interp,
					yerror=σn,
					label="linear interp. obs")	

    Plots.scatter!(t,
					ytrue,
					yerror=σn,
					label="true obs")	


end

# ╔═╡ a9be370b-fe12-40ec-a054-fd864e0694a9
t

# ╔═╡ f513efba-c46d-480d-ac16-13c108aabde6
# use interp instead of expand
# also need the function for control variables (later move into master function)
irregular_sample_control_variables_interp(u) = [irregular_sample_fourier_transform_interp(t[i],FourierTransform(u,ffundamental)) for i in eachindex(t)]

# ╔═╡ 347b5048-68e7-4be6-b386-5b8be36407b3
interp_func = linear_interpolation(collect(xtrue.time), collect(xtrue.x))

# ╔═╡ 346cbf20-ca62-483e-a3eb-905879c38dc6
interp_func.(t)

# ╔═╡ 9e2915e4-ac22-4b27-abf7-755bb1d52ad9
begin   
    # save samples and expected noise together in an `Estimate`
    y_estimate = Estimate(y, fill(σn, M))

    # solve for BLUE of Fourier Transform
    u = BlueFourierTransform.blue_fourier_transform(
        #irregular_sample_control_variables, # takes input and get obs
        irregular_sample_control_variables, # takes input and get obs
        y_estimate, # observations
        Ψ, # first guess frequency spectrum
        1.0) # uncertainty of mean value, not given by spectrum
end

# ╔═╡ 1e5030de-bac9-4679-872a-ca314e2f2f68
beta_tmp = FourierTransform(u.v, ffundamental)

# ╔═╡ 0f91ac16-3bb5-4c53-81a7-3abdcd03724c
xtmp = RegularTimeseries(FourierTransform(u.v, ffundamental))

# ╔═╡ 04c4cb35-b93b-43b1-bf79-2c0130a4d3ae
interp_tmp = linear_interpolation(
            xtmp.time,
            xtmp.x)

# ╔═╡ 06a622cf-483f-428c-a28f-d6b8ce531ae7
interp_tmp(t[1])

# ╔═╡ 5605af27-918a-46bf-bb12-dc8de644958f
collect(xtmp.time)

# ╔═╡ 9718fc83-ba52-467b-bc30-b19e54b2bc6e
irregular_sample_fourier_transform_interp(t, FourierTransform(u.v, ffundamental))

# ╔═╡ da8c32e9-afb0-44b4-bb3a-483f5dc8a51f
irregular_sample_control_variables(u.v)

# ╔═╡ 89437c5a-bfa8-4b17-b1a3-a1ef9136d114
irregular_sample_control_variables_interp(u.v, t)

# ╔═╡ 6a1f9bf7-2f6e-4924-9f8a-265327362668
# how well does it recontruct the obs?
ỹ = irregular_sample_control_variables(u.v)

# ╔═╡ 48149e7c-a344-4654-839d-044fa6c54d3d
function Plots.scatter(xhat::FourierTransform)
	scatter(xhat.freq,
			sqrt.(real.(xhat.coeff).^2 + imag.(xhat.coeff).^2),
			label = "|x̂|",
			markersize = 6)
	scatter!(xhat.freq,
			 real.(xhat.coeff),
			 label = "Real(x̂)", 
		 	legend = :bottomright)
	scatter!(xhat.freq,
			 imag.(xhat.coeff),
			 label = "Imag(x̂)")
	plot!(xlabel = "Dimensional frequency",
		  ylabel = "θ",
		  title="Fourier coefficients")
end

# ╔═╡ 82647927-f8ff-410e-b2dc-9b1290cf2706
scatter(ytrue,ỹ,
	   xlabel="truth",
	   ylabel="obs",
	   label=false)

# ╔═╡ 34e30056-20d4-46c8-87c0-70a21a958e80
Plots.scatter(x̂true)

# ╔═╡ e1c1cc7c-b3db-42fa-b0c4-0f7809ef7437
begin
    # how close is the Fourier Transform description of the data?
    # they actually look pretty different
    x̂̃ = FourierTransform(u.v, ffundamental)
    Δx_real = real.(x̂true.coeff) .- real.(x̂̃.coeff)
    Δx_imag = imag.(x̂true.coeff) .- imag.(x̂̃.coeff)
end

# ╔═╡ 704ac180-7cb5-44e7-a5f1-c0c66463ba11
Plots.scatter(x̂̃)

# ╔═╡ 3ef72948-6090-4d99-9c4d-7c48f2ef4af4
# how faithful is it to the input spectrum?
Ψ̃ = periodogram(x̂̃)    

# ╔═╡ 2891abd0-3f3d-4a02-b70e-756197e4ea3c
Ψ

# ╔═╡ 79cccfe6-334c-4018-a9c2-e294fbab0757
x̂true

# ╔═╡ 294210db-c576-4416-af4b-18e44a5f5c92
x̂̃

# ╔═╡ 612680ed-b595-4137-84a1-f21da7a5eb10
# plot spectra
begin
	
plot(Ψ.freq,
	Ψ.psi,
	#ylims=(1e-5,1e3),
	yscale=:log10,
 	xscale=:log10,
	color=:black,
	lw=4,
	title="frequency spectra",
	label="truth"
	)

	plot!(Ψ̃.freq,
	Ψ̃.psi,
	#ylims=(1e-5,1e3),
	yscale=:log10,
 	xscale=:log10,
	color=:red,
	lw=4,
	label="reconstructed"
	)
end

# ╔═╡ 85916c6d-596a-4fc1-85e2-d6e807298e78
#### Check reconstructed phase
ϕtrue = phase(x̂true)

# ╔═╡ 47cf58a0-f311-4a24-aba1-69ba35a32cd6
ϕ̃ = phase(x̂̃)

# ╔═╡ 3bd45108-1c9c-49b4-933e-c6370136e38f
begin 
	plot(ϕ̃[1:end], label="ϕ̃")
	plot!(ϕtrue[1:end], label="ϕ_true")
	
end

# ╔═╡ Cell order:
# ╠═4b0a896e-b28a-11f1-b73f-77977d857f62
# ╠═8b2fe3cb-4557-4634-82b8-1e14e403bfcb
# ╠═09001ed8-00a9-4d65-bf2f-70a87fac3e4f
# ╠═1e8fe281-32a0-4923-8b1a-b35d6bdc28a7
# ╠═8cbfce62-ea81-42ed-95bd-3bfc052c50ac
# ╠═5ca1234f-374e-4838-834d-96a792125b9a
# ╠═16f76114-2812-4ee0-b6bf-6b9832d08757
# ╠═f7a1395a-b65a-4d84-a062-8851b82af3ad
# ╠═81ea80de-4e29-46e5-a5b2-3ab12faee0b5
# ╠═4d278660-dd7a-4749-9765-67ad369eefa2
# ╠═7ce0b464-a849-4211-8a14-314fb4294eb3
# ╠═6e0e263e-4070-4d00-9769-769c03bc2ef2
# ╠═bcd66b85-4c58-41d5-8958-b79142080f1d
# ╠═8b66a95a-1e69-4ade-9eec-da0a1d7aef2c
# ╠═c746ff66-9c81-48d7-86a8-aee1a2c05dc8
# ╠═dff337bb-9954-4663-8b3f-d34cc16a44b6
# ╠═32724079-1a28-4aa0-8dba-2851471e38ed
# ╠═1e5030de-bac9-4679-872a-ca314e2f2f68
# ╠═0f91ac16-3bb5-4c53-81a7-3abdcd03724c
# ╠═04c4cb35-b93b-43b1-bf79-2c0130a4d3ae
# ╠═5605af27-918a-46bf-bb12-dc8de644958f
# ╠═a9be370b-fe12-40ec-a054-fd864e0694a9
# ╠═06a622cf-483f-428c-a28f-d6b8ce531ae7
# ╠═f513efba-c46d-480d-ac16-13c108aabde6
# ╠═9718fc83-ba52-467b-bc30-b19e54b2bc6e
# ╠═da8c32e9-afb0-44b4-bb3a-483f5dc8a51f
# ╠═89437c5a-bfa8-4b17-b1a3-a1ef9136d114
# ╠═347b5048-68e7-4be6-b386-5b8be36407b3
# ╠═346cbf20-ca62-483e-a3eb-905879c38dc6
# ╠═9e2915e4-ac22-4b27-abf7-755bb1d52ad9
# ╠═6a1f9bf7-2f6e-4924-9f8a-265327362668
# ╠═82647927-f8ff-410e-b2dc-9b1290cf2706
# ╠═34e30056-20d4-46c8-87c0-70a21a958e80
# ╟─48149e7c-a344-4654-839d-044fa6c54d3d
# ╠═704ac180-7cb5-44e7-a5f1-c0c66463ba11
# ╠═e1c1cc7c-b3db-42fa-b0c4-0f7809ef7437
# ╠═3ef72948-6090-4d99-9c4d-7c48f2ef4af4
# ╠═2891abd0-3f3d-4a02-b70e-756197e4ea3c
# ╠═79cccfe6-334c-4018-a9c2-e294fbab0757
# ╠═294210db-c576-4416-af4b-18e44a5f5c92
# ╠═612680ed-b595-4137-84a1-f21da7a5eb10
# ╠═85916c6d-596a-4fc1-85e2-d6e807298e78
# ╠═47cf58a0-f311-4a24-aba1-69ba35a32cd6
# ╠═3bd45108-1c9c-49b4-933e-c6370136e38f
