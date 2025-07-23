include(joinpath(pwd(), "MPC_functions.jl"))
using  Luna
import FFTW                    
import Luna: Hankel  
import NumericalIntegration: integrate, SimpsonEven          
import Dates                   
using DelimitedFiles
using LaTeXStrings
using LsqFit
using Printf
using DSP
using HDF5
using PyPlot

# References
# [1] Marc Hanna, Nonlinear Optics in Multipass Cells, 2021. https://doi.org/10.1002/lpor.202100220
# [2] Anne-Lise Viotti, Multi-pass cells for post-compression of ultrashort laser pulses, 2022, https://doi.org/10.1364/OPTICA.449225 
# [3] https://www.rp-photonics.com
# [4] Saleh, Teich, Fundamentals of Photonics, 2007, Wiley, 2nd Edition, p.54
# [5] https://github.com/LupoLab/Luna.jl/issues/323

# ----------------- SET INPUT FILES AND OUTPUT DIRECTORIES -----------------
    
    # Output directory
    if !isdir(joinpath(@__DIR__, "output"))   # create general output directory if doesnt exist already
        mkdir(joinpath(@__DIR__, "output"))
    end
    out_path = joinpath("output","run_"*Dates.format(Dates.now(), "yyyy_mm_dd__HH_MM_SS"))  # create output subdirectory for this run
    mkdir(out_path)

    # FROG data
    path_FROG_Eω = joinpath("input", "tangerine_FROG_50kHz", "Speck.dat")

    # mirror data
    path_m1 = joinpath("input", "mirrordata_0fs2.txt") # "mirrordata_0fs2.txt"; "mirror_R_100_GDD_0fs2.txt"; "mirror_R_data_GDD_0fs2.txt"
    path_m2 = joinpath("input", "mirrordata_-30fs2.txt") # "mirrordata_-30fs2.txt"; "mirror_R_100_GDD_-30fs2.txt"; "mirror_R_data_GDD_-30fs2.txt"

# ----------------- SET OVERALL SIMULATION OPTIONS -----------------
    use_FROG = false             # whether to use a gaussian or a FROG reconstructed pulse from data as the input    
    use_mirr_data = true        # whether to include the MPC mirror reflectivity and dispersion from given data
    use_taylor_disp = false     # whether to account for the MPC mirror dispersion via a Taylor expansion with coefficients defined below
    
    nl_mode_matching = false    # whether to adjust the MPC modematching to nonlinear focusing
    
    save_intermediate = true    # whether to save plots of intermediate results
    save_gaustic = true         # whether to track the gaustic of each pass
    # save_B_int = false

    comments = " "      # comments to add to the params.txt file for this run

# ----------------- SET PLOTTING OPTIONS -----------------

    res = 300       # resolution (dpi)
    smallval = 1.0e-50      # small value to add in plots where log10 is taken; 
                            # otherwise empty spots appear in plots since PyPlot ignores data points which are 0, since log10(0) = -Inf

# ----------------- SET PHYSICAL PARAMETERS -----------------

    # MPC PARAMETERS
    Nrt = 15                        # number of round trips
    Npass = 2*Nrt                # number of passes
    
    R = 200e-3 #200e-3                      # mirror radius of curvature [m]
    k = Nrt-2                       # for k->Nrt MPC operation is closest to the stability edge and wm is maximized, reducing mirror damage 
    L = R*(1-cos(pi*k/Nrt)) #551e-3 #        # MPC cell length [m] (see Anne-Lise Viotti, Multi-pass cells for post-compression of ultrashort laser pulses, 2022, https://doi.org/10.1364/OPTICA.449225)
    C = L/R

    mtype = :spherical              # mirror type, either :spherical or :parabolic
    ϕm1 = [0,0,0,0]                 # mirror 1 dispersion expressed as a taylor expansion
    ϕm2 = [0,0,0,0]                 # mirror 2 dispersion expressed as a taylor expansion
    
    pres = 1.5 #10e-10 #0.5 #1.5 #10e-10 #1.5    # gas pressure [bar]
    gas = :Kr #:Xe #:Kr                       # gas type
    ion = true                     # if true: enable ionisation response, if false: disable ionisation 
    ion_model="PPT"                 # set to "ADK" or "PPT"; "ADK" is less accurate at low intensities but faster; "PPT" may crash at very high intensities

    # PULSE PARAMETERS
    λ0 = 1030e-9        # central wavelength [m]
    τ = 150e-15 #350e-15 #150e-15         # pulse duration [s]; ignored if FROG spectrum is used
    E_pulse = 250e-6 # 200e-6 #250e-6    # pulse energy [J]
    #w0 = 150e-6         # beam waist at focus [m]      # instead calculated later
    # M2 = 1.0           # beam quality                 # not yet implemented; how would you do that?
    N0, n0, n2 = Tools.getN0n0n2(PhysData.wlfreq(λ0), gas; P=pres, T=PhysData.roomtemp)        # gas number density, linear and nonlinear refractive index [m^2/W]
    
    # calculate mode-matched beam waist for given MPC geometry
    if nl_mode_matching == true
        # see Marc Hanna, "Nonlinear beam matching to gas-filled miltupass cells"

        P_peak = 0.94*E_pulse/τ             # peak power; 0.94 for gaussian pulses (see https://www.rp-photonics.com/peak_power.html)
        P_crit = 3.77*λ0^2/(8*pi*n0*n2)     # critical power; above which beam collapses
        σ = 1-P_peak/(sqrt(2)*P_crit)       # nonlinear correction factor that accounts for self focusing/Kerr lensing

        w0 = sqrt(λ0*L * sqrt(σ*(2*R/L -1)) /(2*pi))   #123.62e-6 #125e-6       # mode-matched waist for nonlinear regime
        wm = w0*sqrt(2*R/(2*R-L)) # not correct I suppose?

    else
        # (see: Anne-Lise Viotti, Multi-pass cells for post-compression of ultrashort laser pulses, 2022, https://doi.org/10.1364/OPTICA.449225 eq. 3&4 
        #  or Marc Hanna, Nonlinear Optics in Multipass Cells, 2021. https://doi.org/10.1002/lpor.202100220 eq. 2&6)
        
        # beam waist in the focus
        w0 = sqrt(R*λ0 * sqrt(C*(2-C)) /(2*pi))         # ! this value is used to define the gaussian beam in the beginning of the simulation
        # w0 = sqrt(λ0*L * sqrt(2*R/L -1) /(2*pi))   #123.62e-6 #125e-6       # these two eq. for w0 should give the same value                                         

        # beam waist on the mirror
        # wm = sqrt(R*λ0 * sqrt(C/(2-C)) /pi) 
        wm = w0*sqrt(2*R/(2*R-L))
    end

    # w0 = 123.57e-6 #123.62e-6 #122.915895e-6
        

# ----------------- SET SIMULATION GRID ----------------------------

    # simulation grid parameters
    λlims = (700e-9, 1400e-9) #(600e-9, 1500e-9) #(700e-9, 1400e-9) # # wavelength range of interest [m]
    trange = 5*τ #2e-12 #5*τ #10*τ #300.0e-15 #0.05e-12    # total extent of time window required [s] (NOTE: if this is too short, the range is extended automatically
    Nz = 201                            # number of points along z at which the spectrum is saved

    # Hankel transformation
    R_Hankel = 3*wm #6e-3 #3.0*wm           # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  
    N_Hankel = 128 #64 #128 #256 #512 #1024             # sample size for Hankel tansform grid 

    q = Hankel.QDHT(R_Hankel, N_Hankel, dim=2)                  # set up discrete Hankel transform matrix, transformation done along 2nd dimension 
    q_1D = Hankel.QDHT(R_Hankel, N_Hankel, dim=1)               # to be applied to arrays that where integrated over ω, so r becomes 1rst dimension
    r = q.r                         # sampled radii [m] from q excluding r=0 so that r = [r1, r2, ..., rn]
    rsym = Hankel.Rsymmetric(q)     # sampled radii [m] mirrored around and including r=0 so that rsym = [–rn, ...-r2, -r1, 0, r1, r2, ..., rn]

    grid = Grid.RealGrid(L, λ0, λlims, trange)               # set up time & space grid for gradient approximation
    ω = grid.ω                      # sampled angular frequencies [rad/s];  NOTE: Nω = (Nt/2) + 1 !
    f = ω./2π                      # sampled linear frequencies [Hz]
    λ = PhysData.wlfreq.(ω)         # sampled wavelengths [m]
    λ[1] = 1.0                        # avoid "Inf" value for DC frequency
    t = grid.t                      # sampled points in time [s];           NOTE: Nt = 2*(Nω-1) !

    energyfun, energyfun_ω = Fields.energyfuncs(grid, q)    # "energyfun" gives total energy in a field E(t); energyfun_ω is needed in case field is defined via FROG spectrum  


# ----------------- SET KERR EFFECT AND/OR PLASMA FORMATION ----------------------------

    ionpot = PhysData.ionisation_potential(gas)                 # set gas ionisation potential   

    if ion_model=="ADK"
        ionrate = Ionisation.ionrate_fun!_ADK(gas)                  # set gas ionisation rate (ADK)
    elseif ion_model=="PPT"
        ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)        # set gas ionisation rate (PPT)
    end 

    # nonlinear response function
    if ion == true  
        # with ionisation and Kerr effect
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                    Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),)  
    elseif ion == false   
        # just Kerr effect
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)    
    end 

# ---------------- PREPARE LUNA SIMULATION -------------------
    if use_FROG==true
        inputs = () # leave input field empty as we will be defined later by the FROG spectrum
    else
        # input field
        # to get a converging beam, we define a Gaussian spot with size w0
        # and then *back-propagate* by L/2 i.e. by half of our propagation window
        # this means the nominal focus is halfway through our propagation in this example
        inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ, energy=E_pulse, w0=w0, propz=-L/2)
    end

    n_gas = PhysData.ref_index_fun(gas, pres)                       # refractive index of the gas for given pressure (no pressure profile, pressure is constant)
    linop = LinearOps.make_const_linop(grid, q, n_gas)              # generate linear operator for pulse-propagation equation
    normfun = NonlinearRHS.const_norm_radial(grid, q, n_gas)        # generate normalisation function for radial symmetry

    densityfun = let dens0=PhysData.density(gas, pres)              # pressure is constant in MPC
        z -> dens0 
    end

    Eωk, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)  # set up simulation for one pass; Eωk is overwritten when FROG spectrum is used


# ---------------- DEFINE MIRRORS ----------------

    # MIRROR 1
        data = readdlm(path_m1, skipstart=1)

        λ_m1 = data[:,1].*1e-9     # read in wavelengths [m]
        R_m1 = data[:,2].*1e-2     # read in reflectivity [fractions]
        GDD_m1 = data[:,3].*1e-30  # read in GDD [s^2]

        # for planck taper function
        windowwidth1 = 20e-9
        λ_min1 = minimum(λ_m1)+windowwidth1
        λ_max1 = maximum(λ_m1)-windowwidth1

        # add additional group delay if necessary
        delay1 = 198.7e-15


    # MIRROR 2
        data = readdlm(path_m2, skipstart=1)

        λ_m2 = data[:,1].*1e-9     # read in wavelengths [m]
        R_m2 = data[:,2].*1e-2     # read in reflectivity [fractions]
        GDD_m2 = data[:,3].*1e-30  # read in GDD [s^2]

        # for planck taper function
        windowwidth2 = 20e-9
        λ_min2 = minimum(λ_m2)+windowwidth2
        λ_max2 = maximum(λ_m2)-windowwidth2

        # add additional group delay if necessary
        delay2 = nothing



# ----------------- DEFINE INPUT VIA COMPLEX FROG SPECTRUM -----------------
    if use_FROG==true

        println("using FROG data as input")

        # READ IN AND PREPARE FROG DATA 
            data = readdlm(path_FROG_Eω)

            λdat = data[:,1]     # read in wavelengths [nm]
            indices = findall(x -> x >= 0, λdat)     # cut off negative wavelengths
            λdat = λdat[indices]
            ωdat = PhysData.wlfreq.(λdat*1e-9)  # calculate angular frequencies 

            E_real = data[indices,4]   # read in real part of complex field
            E_imag = data[indices,5]   # read in imaginary part of complex field
            Eωdat = E_real .+ 1im.*E_imag   # combine real and imaginary part

        # CREATE DATAFIELD (taken from https://github.com/LupoLab/Luna.jl/issues/357)
            df = Fields.DataField(ωdat, Eωdat; energy=1) # energy is irrelevant here as we will rescale again
            Eω = df(grid, nothing) # second argument is unused for DataField (need the planned FT for other types)

            Eωr = Eω .* sqrt.(Maths.gauss.(q.r, w0/2))' # combine with spatial distribution; beam waist taken from above calculation

            Eωk = q * Eωr                              # transform to reciprocal space
            Eωk .*= sqrt(E_pulse)/sqrt(energyfun_ω(Eωk)) # rescale to actual energy

            Fields.prop!(Eωk, -L/2, grid, q)           # from the waist linearly propagate back by -L/2 (to mirror position)
    end

# ----------------- RUN SIMULATION ----------------------------

function run_MPC_simulation(Eωk, output_initial, Npass)
    
    start = Dates.now()
    printstyled("Starting simulation\n", bold=true, color=:green)
    
    output = output_initial
    
    for pass in 1:Npass

        printstyled("pass number: $(pass) / $(Npass)\n", bold=true, color=:green)

        Luna.run(Eωk, grid, linop, transform, FT, output)  # run simulation for one pass 

        zout = output.data["z"]     # z-steps 
        Eωk = output.data["Eω"]     # spectrum for all z-steps

        start_Hankel = Dates.now()
            Eωr = q \ Eωk               # convert to real space
        time_Hankel = Dates.now()-start_Hankel
        tHstring = Utils.format_elapsed(time_Hankel)
        @printf("Hankel transform finished in %s\n", tHstring)

        Eωr_end = Eωr[:, :, end]    # spectrum for last z-step; real space
 
        # track spectral and temporal evolution for on-axis sample (r=0)
        Eωr0_pass = dropdims(Hankel.onaxis(Eωk, q), dims=2)  # spectrum for r=0
        Eωr0[:,:,pass] = Eωr0_pass                    # save output spectrum at r=0 of each pass
 
        if pass==1
            Eωr_in[:, :] = Eωr[:, :, 1]    # input spectrum
        elseif pass==Npass
            Eωr_out[:, :] = Eωr[:, :, end]    # output spectrum
        end

        if save_intermediate==true
            plot_intermediate(Eωr0_pass, pass)
        end
        
        if save_gaustic==true
            
            Iωr = abs2.(Eωr)    # convert to intensity

            # integrate intensity over frequency axis
            Ir = zeros(size(Iωr,2), size(Iωr,3))
            for ii = 1:size(Iωr, 3)
                Ir[:, ii] = integrate(ω, Iωr[:, :, ii], SimpsonEven());
            end

            # create version of Ir mirrored across x-axis
            Ir_sym = Hankel.symmetric(Ir, q_1D)

            # Calculate beam waist for each z-step
            w_z_pass = vec(Maths.fwhm(rsym, Ir_sym, level=1/exp(2), method=:spline) /2)
            append!(w_z, w_z_pass) #caution: appended version will have a double entry for the first and last beam waist of each pass

            # Calculate minimum beam waist and beam waist at mirror position for each pass
            push!(w_min, minimum(w_z_pass))
            push!(w_mirr, w_z_pass[end])
        end

        refocus!(Eωr_end; mtype=mtype)  # refocus diverging beam
        
        # apply mirror dispersion and reflectivity from data alternately
        if use_mirr_data==true
            if isodd(pass)
                reflecting_mirror1!(Eωr_end) 
            else
                reflecting_mirror2!(Eωr_end) 
            end
        end

        # apply mirror dispersion directly via Taylorcoeff. of the spectral phase
        if use_taylor_disp==true
            if iseven(pass)
                Fields.prop_taylor!(Eωr_end, grid, ϕm1, λ0)
            else
                Fields.prop_taylor!(Eωr_end, grid, ϕm2, λ0)
            end
        end

        Eωk = q*Eωr_end  # transform back to k-space

        # initialize new MemoryOutput except for last pass
        if pass==Npass
            break
        else
            output = Output.MemoryOutput(0, grid.zmax, Nz)
        end

    end

    totaltime = Dates.now()-start
    dtstring = Utils.format_elapsed(totaltime)
    printstyled(@sprintf("Simulation finished in %s\n", dtstring); bold=true, color=:green)

    return output
end

# preallocate arrays to store results
Eωr_in  = Array{ComplexF64,2}(undef, length(ω), length(r))  # input spectrum
Eωr_out = Array{ComplexF64,2}(undef, length(ω), length(r))  # output spectrum
Eωr0    = Array{ComplexF64,3}(undef, length(ω), Nz, Npass)  # spectrum for all z-steps of all passes

w_z = Float64[]        # beam waist for all z-steps over all passes
w_min = Float64[]      # minimum beam waist per pass (focus spot size)    
w_mirr = Float64[]     # beam waist at the end of each pass (mirror spot size)

println("R = ", R)
println("L = ", L)
println("L/R = ", C)
println("n2 = ", n2)
# println("σ = ",σ)
@printf("Stationary focal beam waist w0 = %.2f µm\n", w0*1e6)
@printf("Stationary mirror beam waist wm = %.2f µm\n", wm*1e6)

output_initial = Output.MemoryOutput(0, grid.zmax, Nz)    # configure initial output
output = run_MPC_simulation(Eωk, output_initial, Npass)   # run simulation


# ----------------- SAVE RESULTS ----------------------------
print("Saving data...")

z = output.data["z"]     # z-steps
Eωk = output.data["Eω"]

open(joinpath(out_path,"params.txt"), "w") do file

    write(file, "# MPC parameters\n")
    write(file, "L          = $(L)\n")
    write(file, "R          = $(R)\n")
    write(file, "L/R        = $(C)\n")
    write(file, "k          = $(k)\n")
    #write(file, "Nrt        = $(Nrt)\n")       # disabled for now, since I am using the number of passes for testing
    write(file, "Npass      = $(Npass)\n")
    write(file, "mtype      = $(mtype)\n")
    write(file, "propz      = $(-L/2)\n")
    write(file, "\n")

    write(file, "gas        = $(gas)\n")
    write(file, "pres       = $(pres)\n")
    write(file, "n0         = $(n0)\n")
    write(file, "n2         = $(n2)\n")
    write(file, "ion        = $(ion)\n")
    if ion write(file, "ion_model  = $(ion_model)\n") end
    write(file, "\n")


    write(file, "# pulse parameters\n")
    write(file, "λ0         = $(λ0)\n")
    if use_FROG==false 
        write(file, "τ          = $(τ)\n") 
    end
    write(file, "E_pulse    = $(E_pulse)\n")
    write(file, "w0         = $(w0)\n")
    write(file, "wm         = $(wm)\n")
    # write(file, "M2         = $(M2)\n")  
    # write(file, "σ          = $(σ)\n")
    if use_FROG == true
        write(file, "pulse shape    = from FROG reconstruction\n") 
        write(file, "path_FROG_Eω   = $(path_FROG_Eω)\n")
    else
        write(file, "pulse shape    = gaussian\n") 
    end      
    write(file, "\n")


    write(file, "# simulation grid parameters\n")
    write(file, "λlims      = $(λlims)\n")
    write(file, "trange request = $(trange)\n")
    write(file, "trange used    = $(t[end]-t[1])\n")
    write(file, "Nω         = $(length(ω))\n")
    write(file, "Nt         = $(length(t))\n")
    write(file, "Nz         = $(Nz)\n")
    write(file, "\n")

    write(file, "# Hankel transformation\n")
    write(file, "R_Hankel   = $(R_Hankel)\n")
    write(file, "N_Hankel   = $(N_Hankel)\n")
    write(file, "\n")

    write(file, "# MPC mirrors\n")
    if use_mirr_data==true
    write(file, "path_mirror1 = $(path_m1)\n")
    write(file, "windowwidth1 = $(windowwidth1)\n")
    write(file, "λ_min1     = $(λ_min1)\n")
    write(file, "λ_max1     = $(λ_max1)\n")
    write(file, "delay1     = $(delay1)\n")
    write(file, "path_mirror2 = $(path_m2)\n")
    write(file, "windowwidth2 = $(windowwidth2)\n")
    write(file, "λ_min2     = $(λ_min2)\n")
    write(file, "λ_max2     = $(λ_max2)\n")
    write(file, "delay2     = $(delay2)\n")
    else
    write(file, "no mirror dispersion or reflectivity used from data.\n")
    end

    if use_taylor_disp==true
    write(file, "ϕm1     = $(ϕm1)\n")
    write(file, "ϕm2     = $(ϕm2)\n")
    else
    write(file, "no taylor expansion used for mirror dispersion.\n")
    end
    write(file, "\n")

    
    write(file, "comments   = $(comments)\n")
end

params = Dict(
    "λ0" => λ0,
    "τ" => τ,
    "E_pulse" => E_pulse,
    "w0" => w0,
    "wm" => wm,
    # "M2" => M2,
    # "σ" => σ,
    "L" => L,
    "R" => R,
    "C" => C,
    "k" => k,
    "Npass" => Npass,
    "mtype" => mtype,
    "propz" => -L/2,
    "gas" => gas,
    "pres" => pres,
    "n0" => n0,
    "n2" => n2,
    "ion" => ion,
    "ion_model" => ion_model,
    "λmin" => λlims[1],
    "λmax" => λlims[2],
    "trange_request" => trange,
    "trange_used" => t[end]-t[1],
    "Nω" => length(ω),
    "Nt" => length(t),
    "Nz" => Nz,
    "R_Hankel" => R_Hankel,
    "N_Hankel" => N_Hankel
)


h5open(joinpath(out_path, "output.h5"), "w") do file
    
    # Save main data; complex electric field in frequency domain of each pass
    dset = create_dataset(
        file, 
        "Eωr0", 
        datatype(Eωr0), 
        size(Eωr0); 
        chunk=(64,32,1), # tune as needed
        compress=6
    )
    write(dset, Eωr0)

    # Save in and output spectrum (with r dependecy)
    file["Eωr_in"] = Eωr_in    # input spectrum
    file["Eωr_out"] = Eωr_out    # output spectrum

    # # Save beam waists
    file["w_z"] = w_z
    file["w_min"] = w_min
    file["w_mirr"] = w_mirr

    # Save axes 
    file["ω"] = ω
    file["t"] = t
    file["r"] = r
    file["z"] = z

    # Save parameters in a group
    param_group = create_group(file, "params")
    for (key, val) in params
        # HDF5 does not support writing general objects as attributes.
        # Save as scalar or 1D arrays if supported, else convert to string
        try
            param_group[string(key)] = val
        catch e
            # fallback for e.g. NamedTuple, custom types etc.
            param_group[string(key)] = string(val)
        end
    end
end

println("done.")
