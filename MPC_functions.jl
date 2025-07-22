function plot_intermediate(Eωr0_pass, pass)
            
    ϕω_in = unwrap_phase(ω, t, Eωr0_pass[:,1])
    ϕω_out = unwrap_phase(ω, t, Eωr0_pass[:,end])

    ϕω_in = blank_phase(ω, abs2.(Eωr0_pass[:,1]), ϕω_in; level=0.05)
    ϕω_out = blank_phase(ω, abs2.(Eωr0_pass[:,end]), ϕω_out; level=0.05)

    # Input spectrum
    fig, ax1 = plt.subplots()
    ax1.set_title("Pass $pass: input spectrum")
    ax1.plot(λ.*1e9, Maths.normbymax(abs2.(Eωr0_pass[:,1])), label="Intensity")
    ax2 = ax1.twinx()
    ax2.plot(λ.*1e9, ϕω_in, label="ϕω", color="orange")
    ax1.set_xlim(λlims[1]*1e9, λlims[2]*1e9)
    ax1.set_xlabel("λ [nm]")
    ax1.set_ylabel("I [norm.]")
    ax2.set_ylabel("phase [rad]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    plt.tight_layout()     
    plt.savefig(joinpath(out_path,"spectrum_pass_$(pass)_input.png"),dpi=res)
    plt.close()

    # Output spectrum
    fig, ax1 = plt.subplots()
    ax1.set_title("Pass $pass: output spectrum")
    ax1.plot(λ.*1e9, Maths.normbymax(abs2.(Eωr0_pass[:,end])), label="Intensity")
    ax2 = ax1.twinx()
    ax2.plot(λ.*1e9, ϕω_out, label="ϕω", color="orange")
    ax1.set_xlim(λlims[1]*1e9, λlims[2]*1e9)
    ax1.set_xlabel("λ [nm]")
    ax1.set_ylabel("I [norm.]")
    ax2.set_ylabel("phase [rad]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    plt.tight_layout()     
    plt.savefig(joinpath(out_path,"spectrum_pass_$(pass)_output.png"),dpi=res)
    plt.close()


    Etr0_in = Maths.hilbert(FFTW.irfft(Eωr0_pass[:,1], length(t), 1))  
    Etr0_out = Maths.hilbert(FFTW.irfft(Eωr0_pass[:,end], length(t), 1))  

    ϕt_in = unwrap(angle.(Etr0_in))
    ϕt_out = unwrap(angle.(Etr0_out))

    ϕt_in = blank_phase(t, abs2.(Etr0_in), ϕt_in; level=0.05)
    ϕt_out = blank_phase(t, abs2.(Etr0_out), ϕt_out; level=0.05)

    # Input pulse
    fig, ax1 = plt.subplots()
    ax1.set_title("Pass $pass: input pulse")
    ax1.plot(t.*1e15, Maths.normbymax(abs2.(Etr0_in)), label="Intensity")
    ax2 = ax1.twinx()
    ax2.plot(t.*1e15, detrend(ϕt_in), label="ϕt", color="orange")
    ax1.set_xlabel("t [fs]")
    ax1.set_ylabel("I [norm.]")
    ax2.set_ylabel("phase [rad]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    ax1.set_xlim(minimum(t)*1e15, maximum(t)*1e15)
    plt.tight_layout()
    plt.gcf()     
    plt.savefig(joinpath(out_path,"pulse_pass_$(pass)_input.png"),dpi=res)
    plt.close()

    # Output pulse
    fig, ax1 = plt.subplots()
    ax1.set_title("Pass $pass: output pulse")
    ax1.plot(t.*1e15, Maths.normbymax(abs2.(Etr0_out)), label="Intensity")
    ax2 = ax1.twinx()
    ax2.plot(t.*1e15, detrend(ϕt_out), label="ϕt", color="orange")
    ax1.set_xlabel("t [fs]")
    ax1.set_ylabel("I [norm.]")
    ax2.set_ylabel("phase [rad]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    ax1.set_xlim(minimum(t)*1e15, maximum(t)*1e15)
    plt.tight_layout()
    plt.gcf()    
    plt.savefig(joinpath(out_path,"pulse_pass_$(pass)_output.png"),dpi=res)
    plt.close()
end

"""
Removes linear slope. CHECK AGAIN!
"""
# function detrend(x::AbstractVector)
#     n = length(x)
#     t = 1:n
#     X = [ones(n) t]                # design matrix for linear regression: y = a + b*t
#     β = X \ x                      # Least squares solution
#     trend = X * β
#     return x .- trend
# end

function detrend(x::AbstractVector)
    isgood = .!isnan.(x)       # remove possible NaN values from phase_blanking()
    n = length(x)
    t = 1:n

    x_good = x[isgood]
    t_good = t[isgood]

    X = [ones(length(x_good)) t_good]        # Design matrix for non-NaN positions
    β = X \ x_good                          # Fit only to non-NaN values
    trend_good = X * β                      # Trend for non-NaNs

    detrended = similar(x)
    detrended .= NaN
    detrended[isgood] = x_good .- trend_good
    return detrended
end

"""
Change from Iω to Iλ or the other way around.
"""
function Iwlfreq(ωλ, Iωλ)
    λω = PhysData.wlfreq(ωλ)
    Iλω = @. Iωλ*(2*pi*PhysData.c)/(λω^2)
    return Iλω
end

function gauss_fit(t, intensity)
    
    ### normalize intensity (otherwise fit is way off)
    intensity_norm = intensity ./ maximum(intensity)


    ### Define gaussian (fit model)
    function gauss(t, p)
        A, µ, fwhm = p

        σ = @. fwhm / (2*sqrt(2*log(2)))
        I = @. A * exp(-(t - µ)^2 / (2*σ^2))
        
        return I
    end

    ### setting initial guess values for gauss parameters
    A_init = maximum(intensity_norm)                            # amplitude is the maximum intensity
    µ_init = t[argmax(intensity_norm)]                          # shift in x direction is the time with index of maximum intensity
    fwhm_init = Maths.fwhm(t, intensity_norm, method=:nearest, minmax=:max)      # using luna guess for fwhm

    guess = [A_init, µ_init, fwhm_init]

    # performing gaussian fit
    fit = curve_fit(gauss, t, intensity_norm, guess)

    # extracting fitted parameters
    fitted_params = fit.param
    A, µ, fwhm = fitted_params

    # calculate fitted intensities (converting back to input magnitude)
    fitted_intensities = gauss(t, fitted_params) .* maximum(intensity)

    return fwhm, fitted_intensities
end

"""
    unwrap_phase(E_complex; linear=false)

(the following is based on the documentation of the overlap() function in the Modes.jl module of luna)

Pulses centred on t=0 have large linear spectral phase components which can confuse
phase unwrapping and lead to oscillations in the spectral phase.
Here we shift the pulse to t=-t_max of the time grid to remove this, then do the
unwrapping.
The unwrapped phase is returned by default without the linear contribution.
The linear contribution can be put back on the phase by setting "linear=true".

'E_complex' is the complex electric field from which to derive the phase values
'linear'    whether to reintroduce the eliminated large linear contribution to the phase
"""
function unwrap_phase(ω, t, E_complex; linear=false)

    E_unwrap = copy(E_complex)
    
    Nt = length(t)
    dt = t[2] - t[1]

    τ_shift = -Nt*dt/2
    E_unwrap .*= exp.(-1im .*ω .*τ_shift)      # eliminating large linear contribution by shifting to -t_max before unwrapping

    ϕω = unwrap(-angle.(E_unwrap); dims=1)     

    if linear == true
        ϕω .-= ω .*τ_shift      # putting back linear part of the phase
    end

    return ϕω
end



"""
    blank_phase(x, I, ϕ; level=0.05)

Returns a blanked version of the phase array ϕ, meaning all phase values are put to NaN for which the intensity I is below a certain threshold.
'x' is the axis, either time or frequency axis.
'I' are the intensities.
'ϕ' are the phase values.
'level' defines the threshold as a fraction of the maximum intensity above which phase values are not blanked.
"""
function blank_phase(x, I, ϕ; level=0.05)
    
    # works for 1D and 2D I, ϕ arrays
    
    if size(I) != size(ϕ)
        error("Error: I and ϕ have to be of the same size.")
    end
    
    num_col = size(I, 2)  # Number of columns in I and ϕ
    ϕ_blanked = similar(ϕ)  # Allocate array to store results

    # Iterate over each column in I
    for i in 1:num_col
        I_col = I[:, i]
        ϕ_col = ϕ[:, i]

        val = level * maximum(I_col)    # threshold intensity
        maxidx = argmax(I_col)
        xmax = x[maxidx]               

        # find left and right indices between which intensity is above
        # lefti = findlast((x .< xmax) .& (I_col .< val))
        # righti = findfirst((x .> xmax) .& (I_col .< val))

        lefti = findfirst((x .< xmax) .& (I_col .> val))
        righti = findlast((x .> xmax) .& (I_col .> val))

        ϕ_blanked_col = copy(ϕ_col)
        
        if !isnothing(lefti)
            ϕ_blanked_col[1:lefti-1] .= NaN
        end
        if !isnothing(righti)
            ϕ_blanked_col[righti+1:end] .= NaN
        end

        # shift phase values to zero (convinient for comparing the phase graphs)
        ϕ_blanked_col .-= minimum(filter(!isnan, ϕ_blanked_col))

        ϕ_blanked[:, i] = ϕ_blanked_col
    end
    
    return ϕ_blanked
end


function refocus!(Eωr; mtype=:spherical)
    # Refocus the beam by multiplying the electric field with the right phase factor
    # See https://github.com/LupoLab/Luna.jl/issues/323
    
    print("Refocusing beam...")

    # calculate kz(ω, k_perp)
    kzsq = (grid.ω/PhysData.c).^2 .- reshape(q.k.^2, (1, length(q.k)))
    kzsq[kzsq .< 0] .= 0
    kz = sqrt.(kzsq)

    # R is the focal length of the mirror
    # q.r is the radial coordinate in the simulation

    if mtype == :parabolic
        # phase imparted by reflection off of a parabolic mirror at normal incidence
        φr_foc = -kz .* reshape(q.r.^2/R, (1, length(q.r))) 

    elseif mtype == :spherical
        # phase imparted by reflection off of a spherical mirror at normal incidence 
        # (see: Saleh, Teich, Fundamentals of Photonics, 2007, Wiley, 2nd Edition, p.54)
        φr_foc = -kz .* 2 .* reshape(R .- sqrt.(R^2 .- (q.r .^2)), (1, length(q.r)))

    else
        error("Error: mirror type has to be either :spherical or :parabolic")
    end

    #multiply electric field by phase factor
    Eωr .*= exp.(-1im .* φr_foc)    # the '.' ensures that the input Eωr itself is modified and not returned as a new array

    println("done.")

end


function reflecting_mirror1!(Eωr)

    print("Applying reflectivity and GDD data of mirror 1...")
    Fields.prop_mirror!(Eωr, grid.ω, 1, λ_m1, R_m1, λ_m1, GDD_m1, λ0, λ_min1, λ_max1, windowwidth=windowwidth1) 
    
    if delay1 !== nothing
        # adding additional linear phase
        ω0 = PhysData.wlfreq(λ0)
        ϕ_delay = delay1 .* (grid.ω .- ω0)
        Eωr .*= exp.(-1im*ϕ_delay)
    end

    println("done.")
end

function reflecting_mirror2!(Eωr)

    print("Applying reflectivity and GDD data of mirror 2...")
    Fields.prop_mirror!(Eωr, grid.ω, 1, λ_m2, R_m2, λ_m2, GDD_m2, λ0, λ_min2, λ_max2, windowwidth=windowwidth2) 

    if delay2 !== nothing
        # adding additional linear phase
        ω0 = PhysData.wlfreq(λ0)
        ϕ_delay = delay2 .* (grid.ω .- ω0)
        Eωr .*= exp.(-1im*ϕ_delay)
    end

    println("done.")
end