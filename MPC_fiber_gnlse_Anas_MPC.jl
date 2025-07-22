using Luna
using LaTeXStrings
using FFTW
using DelimitedFiles
using FFTW
using DSP

@eval import PyPlot: pygui, plt, PyDict, matplotlib 
    close("all")
    pygui(true)


function ManualPulse(ω, t, Eω; energy=E_pulse, scale_energy=nothing, kwargs...)

    τ = length(t) * (t[2] - t[1])/2 # middle of old time window

    Eωm = Eω[:, end]
    return DataPulse(ω, Eωm .* exp.(1im .* ω .* τ); energy=energy, kwargs...)
end


function find_closest_index(array, value)
    differences = abs.(array .- value)
    return argmin(differences)
end

function plot_spectrum(ω, Eω; phase_range=(minimum(ω), maximum(ω)), components=[:abs2, :unwrapped_phase])
    
    fig, ax1 = plt.subplots()
    ax1.set_title("Complex spectrum")

    idx_min = find_closest_index(ω, phase_range[1])
    idx_max = find_closest_index(ω, phase_range[2])


    # Plot within the first axes based on specified components
    if :abs2 in components
        ax1.plot(ω .* 1e-15, abs2.(Eω), label="abs2", marker="o")
    end

    if :real in components
        ax1.plot(ω .* 1e-15, real.(Eω), label="real", marker="o")
    end

    if :imaginary in components
        ax1.plot(ω .* 1e-15, imag.(Eω), label="imaginary", marker="o")
    end

    # Create a secondary y-axis for phase-related plots
    ax2 = ax1.twinx()

    if :phase in components
        ax2.plot(ω[idx_min:idx_max] .* 1e-15, angle.(Eω[idx_min:idx_max]), label="phase", color="orange")
    end

    if :unwrapped_phase in components
        ax2.plot(ω[idx_min:idx_max] .* 1e-15, unwrap(angle.(Eω[idx_min:idx_max])), label="phase unwrapped", color="orange")
    end

    # Configure axes and layout
    ax1.set_xlabel("ω [PHz]")
    ax1.set_ylabel("a.u.")
    ax2.set_ylabel("phase [rad]")

    # Set legends for both axes
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")

    plt.tight_layout()
    
    return fig, ax1, ax2
end

function mirror!(GDD="low")
    if GDD=="low"
        Fields.prop_mirror!(Eω, ω, 1, λ_low_GDD, R_low_GDD, λ_low_GDD, GDD_low_GDD, λ0, 400e-9, 1800e-9)
    elseif GDD=="high"
        Fields.prop_mirror!(Eω, ω, 1, λ_high_GDD, R_high_GDD, λ_high_GDD, GDD_high_GDD, λ0, 400e-9, 1800e-9)
    else
        println("GDD has to be either \"low\" or \"high\"")
    end
end

function run_split_simulation(datapulse, Npass, γ, L, βs, λ0, λlims, trange)

    output = nothing 

    for pass in 1:Npass
        println("pass number: ", pass)

        output = prop_gnlse(γ, L, βs; λ0, λlims, trange, pulses=datapulse)
        println("next simulation step was successfull.")

        # Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=λlims, trange=(-10*τ, 10*τ))
        # plt.title("Pass number :$(pass)")
        # plt.tight_layout()
    
        ω = output["grid"]["ω"]
        ω = fftshift(ω)
        println("size(ω) = ",size(ω))

        #indices = findall(x -> x >= 0, ω)     # omit negative wavelengths
        #ω = ω[indices]
        #println("size(ω) = ",size(ω))

        Eω = output["Eω"][:,end]
        Eω = fftshift(Eω)
        #Eω = Eω[indices]
        println("size(Eω) = ",size(Eω))

        t = output["grid"]["t"]
        println("size(t) = ",size(t))

        # plt.figure("time axis")
        # plt.title("time axis (unshifted)")
        # plt.plot(t.*1e15)
        # plt.xlabel("index")
        # plt.ylabel("t [fs]")
        # plt.show()

        tshift = length(t) * (t[2] - t[1])/2 # middle of old time window
        Eω = Eω .* exp.(1im .* ω .* tshift)

        # fig, ax1, ax2 = plot_spectrum(ω, Eω, phase_range=(1.78e15, 1.88e15))
        # ax1.set_title("spectrum after pass number $(pass)")
        # #ax1.set_xlim(1.78, 1.88)

        # zslice = [0,L]
        # Plotting.spec_1D(output, zslice; λrange=(850e-9, 1200e-9), log10=false)
        # Plotting.add_fwhm_legends(plt.gca(), "nm")
        # plt.title("M2(input) = $(ceil(M2 * 100) / 100)")
        # plt.tight_layout()


        z = output["z"]
        println("size(z) = ",size(z))

        energy = Processing.energy(output; bandpass=nothing)
        println("size(energy) = ",size(energy))
        println("energy[end] = ", energy[end])


        # apply refectivity and GDD of MPC mirrors
        if use_mirrors==true
            if iseven(pass)
                Fields.prop_mirror!(Eω, ω, 1, λ_low_GDD, R_low_GDD, λ_low_GDD, GDD_low_GDD, λ0, λwindow[1], λwindow[2]) #921e-9, 1168e-9) #800e-9, 1300e-9)#1000e-9, 1050e-9)#)
                
                # Compensate linear time shift
                ω0 = PhysData.wlfreq(λ0)
                ϕ_delay = 200e-15 .* (ω .- ω0)
                Eω .*= exp.(-1im*ϕ_delay)
            else
                Fields.prop_mirror!(Eω, ω, 1, λ_high_GDD, R_high_GDD, λ_high_GDD, GDD_high_GDD, λ0, λwindow[1], λwindow[2]) #921e-9, 1168e-9) #800e-9, 1300e-9)#1000e-9, 1050e-9)#921e-9, 1168e-9)
            end
        end

        datapulse = Pulses.DataPulse(ω, Eω; mode=:lowest, polarisation=:linear, propagator=nothing, energy=energy[end])
        println("creating new datapulse from previous simulation was successfull.")
    end

    return output
end


# SIMULATION OPTIONS
    use_FROG = false
    use_gauss = true

    use_split_sim = true
    use_mirrors = true



# DATA FILES

    # FROG data
    #dir = "tangerine_FROG_1kHz"
    dir = "tangerine_FROG_50kHz"

    file_FROG_Et = "Ek.dat"

    file_FROG_Eω = "Speck.dat"
    path_FROG_Eω = joinpath(dir, file_FROG_Eω)

    # file_FROG_spec = "Speck.dat"                   # name of complex IR FROG input spectrum file 
    # file_FROG_pulse = "20220505_retrievedpulse.txt"


    # mirror data
    file_low_GDD = joinpath("input", "mirrordata_0fs2.txt")
    file_high_GDD = joinpath("input" , "mirrordata_-30fs2.txt")


    # measured MPC spectrum
    file_MPC_spectrum = joinpath("input", "MPC_measured_spectra_2024-10-25", "Spectrum", "Spectrum before fibre at 50kHz.txt")

    num_header_lines = 14
    data = readdlm(file_MPC_spectrum, '\t', Float64, '\n', skipstart=num_header_lines)

    λ_exp = data[:, 1]
    I_exp = data[:, 2]


# SETTING PARAMETERS

    # pulse parameters
    λ0 = 1030e-9        # central wavelength in nm
    τ = 150e-15         # pulse duration 
    E_pulse = 250e-6    # pulse energy 
    M2 = 1.16           # beam quality

    propz = 3e-3 #?  propagation distance from the waist?

    # MPC parameters
    L = 380.19e-3          # MPC cell length in m
    R = 200e-3          # radius of curvature in m
    Nrt = 15           # number of round trips
    Npass =  2*Nrt       # number of passes

    pres = 1.5               # gas pressure in bar
    gas = :Kr           # gas type


    # simulation grid parameters
    λlims = (400e-9, 1800e-9) #(800e-9, 1200e-9) # # wavelength range of interest (?) in nm
    
    #N_grid = 2048                   # sampling size of time axis (?)
    trange = 10*τ #N_grid*dt         # time window 
    #dt = trange/N_grid #1.0e-15     # time step 


    # window for applying mirror data
    λwindow = (921e-9, 1168e-9) #(1000e-9, 1060e-9) #(921e-9, 1168e-9) #(800e-9, 1200e-9)
    


# CALCULATE INPUT FOR LUNA FUNCTION

    flength = 2*Nrt*L     # equivalent waveguide length

    ω0 = PhysData.wlfreq(λ0)    # carrier frequency
    n2 = Tools.getN0n0n2(ω0, gas; P=pres, T=PhysData.roomtemp)[3]        # nonlinear index in [m^2/W]
    println("n2 = ", n2)

    A_eff = M2*λ0*L / (2*atan(sqrt(L/(2R-L))))     # effective area of equivalent waveguide in [m^2]
    println("A_eff = ",A_eff)

    k0 = 2π/λ0          # wavenumber of central wavelength [rad/m]
    γ = k0*n2/A_eff     # nonlinear coefficient in [rad/(m*W)]
    println("γ = ",γ)


# Calculate β for given gas & pressure
    orders = 0:1:7
    βs = zeros(length(orders)) 

    for (i, order) in enumerate(orders)
        βs[i] = PhysData.dispersion(order, gas, λ0, pres) 
    end
    println("βs[3] luna: ",βs[3])


    # in case mirror data is not used: apply mirror GDD directly to overall dispersion of equivalent waveguide
    if use_mirrors==false
        GDD = -30.0e-30/2       # averaged GDD acquired per pass due to chirped Mirrors
        println("GDD considered per pass = ",GDD/L)
        βs[3] += GDD/L
        println("βs[3] chirped: ",βs[3])
    end

# ADD MIRRORS 
if use_mirrors==true

    # LOW GDD MIRROR
    data = readdlm(file_low_GDD, skipstart=1)

    λ_low_GDD = data[:,1].*1e-9     # read in wavelength [m]

    R_low_GDD = data[:,2].*1e-2     # read in reflectivity [fractions]
    GDD_low_GDD = data[:,3].*1e-30  # read in GDD [s]

    fig, ax1 = plt.subplots()
    ax1.set_title("Low GDD mirror")
    ax1.plot(λ_low_GDD.*1e9, R_low_GDD.*1e2, label="reflectivity")#, marker="o")
    ax2 = ax1.twinx()
    ax2.plot(λ_low_GDD.*1e9, GDD_low_GDD.*1e30, label="GDD", color="orange")
    ax2.axhline(y=0.0, color="black", linestyle="--")
    ax2.set_ylim(-50, 50)
    ax1.set_xlabel("λ [nm]")
    ax1.set_ylabel("reflectivity [%]")
    ax2.set_ylabel("GDD [fs2]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    plt.tight_layout()
    plt.show()

    
    # HIGH GDD MIRROR
    data = readdlm(file_high_GDD, skipstart=1)

    λ_high_GDD = data[:,1].*1e-9    # read in wavelength [m]

    R_high_GDD = data[:,2].*1e-2    # read in reflectivity [fractions]   
    GDD_high_GDD = data[:,3].*1e-30 # read in GDD [s]

    fig, ax1 = plt.subplots()
    ax1.set_title("High GDD mirror")
    ax1.plot(λ_high_GDD.*1e9, R_high_GDD.*1e2, label="reflectivity")#, marker="o")
    ax2 = ax1.twinx()
    ax2.plot(λ_high_GDD.*1e9, GDD_high_GDD.*1e30, label="GDD", color="orange")
    ax2.axhline(y=-30.0, color="black", linestyle="--")
    ax2.set_ylim(-50, 50)
    ax1.set_xlabel("λ [nm]")
    ax1.set_ylabel("reflectivity [%]")
    ax2.set_ylabel("GDD [fs2]")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    plt.tight_layout()
    plt.show()
end


# DEFINE INPUT VIA COMPLEX FROG SPECTRUM
if use_FROG==true

    println("using FROG spectrum as input")

    data = readdlm(path_FROG_Eω)

    # READ IN WAVELENGTHS AND CUT OFF NEGATIVE WAVELENGTHS
    λ = data[:,1]     # read in wavelength

    indices = findall(x -> x >= 0, λ)     # omit negative wavelengths
    λ = λ[indices]
    println("size(λ): ",size(λ))

    # CALCULATE FREQUENCY AXIS
    f = @. PhysData.c/(λ*1e-9)
    ω = 2*pi.*f

    # READ IN COMPLEX FROG SPECTRUM
    E_real = data[indices,4]   # read in real part of complex field
    E_imag = data[indices,5]   # read in imaginary part of complex field

    Eω = E_real .+ 1im.*E_imag
    println("size(Eω_tangerine_FROG): ", size(Eω))

    datapulse = Pulses.DataPulse(ω, Eω; mode=:lowest, polarisation=:linear, propagator=nothing, energy=E_pulse)


    # PLOT INPUT SPECTRUM

    fig, ax1, ax2 = plot_spectrum(ω, Eω; phase_range=(1.78e15, 1.88e15))
    ax1.set_title("Tangerine spectrum (first overall input)")


elseif use_gauss==true
    datapulse = Pulses.GaussPulse(;λ0=λ0, τfwhm=τ, energy=E_pulse)

else
    println("Either use_gauss or use_FROG has to be true.") 
end


# RUN SIMULATION

if use_split_sim==false 

    if use_gauss==true
        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)
    elseif use_FROG==true
        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, pulses=datapulse)
    else
        println("Either use_gauss or use_FROG has to be true.")
    end

elseif use_split_sim==true
        # datapulse used in the function is defined already for gauss or FROG case
        output = run_split_simulation(datapulse, Npass, γ, L, βs, λ0, λlims, trange)
end

    

# PLOT RESULTS

    # specx, Iω = Processing.getIω(output, :λ; specrange=λlims)
    # println(size(specx))
    # println(size(Iω))

    # plt.figure()
    # plt.plot(specx*1e9, Iω[:,end])
    # plt.show()

    # open(joinpath("spectrum_split_sim_FROG.txt"), "w") do file
    #     writedlm(file, zip(specx*1e9, Iω[:,end]))
    # end


    Plotting.stats(output)

    Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(850e-9, 1200e-9), trange=(-trange/2, +trange/2)) #(600e-9, 1500e-9)

    if use_split_sim
        zslice = [L]
    else
        zslice = [flength]
    end
    
    Plotting.spec_1D(output, zslice; λrange=(850e-9, 1200e-9), log10=false, label="simulated") #(600e-9, 1500e-9)
    #Plotting.add_fwhm_legends(plt.gca(), "nm")
    plt.title("MPC output spectrum")
    plt.plot(λ_exp, I_exp, label="measured")
    plt.legend()
    plt.tight_layout()


    Plotting.time_1D(output, zslice; trange=(-trange/2, +trange/2))#, bandpass=(180e-9, 220e-9))
    Plotting.add_fwhm_legends(plt.gca(), "fs")