using  Luna
import Luna.PhysData: wlfreq    
import FFTW                    
import Luna: Hankel  
import NumericalIntegration: integrate, SimpsonEven          
import Dates                   
using  DelimitedFiles
using  LaTeXStrings
using LsqFit
using Printf
using DSP

### Aim of this script is to illustrate how the mirrors affect the pulse. ###


function reflecting_mirror1!(Eωr)

    print("Applying reflectivity and GDD data of mirror 1...")
    Fields.prop_mirror!(Eωr, grid.ω, N_reflections, λ_m1, R_m1, λ_m1, GDD_m1, λ0, 820e-9, 1280e-9) #λlims[1], λlims[2])#1000e-9, 1050e-9)#921e-9, 1168e-9)
    println("done.")

    # adding additional linear phase
    ω0 = PhysData.wlfreq(λ0)
    ϕ_delay = 200e-15 .* (grid.ω .- ω0)
    Eωr .*= exp.(-1im*ϕ_delay * N_reflections)

end

function reflecting_mirror2!(Eωr)

    print("Applying reflectivity and GDD data of mirror 2...")
    Fields.prop_mirror!(Eωr, grid.ω, N_reflections, λ_m2, R_m2, λ_m2, GDD_m2, λ0, 820e-9, 1280e-9) #λlims[1], λlims[2])#1000e-9, 1050e-9)#921e-9, 1168e-9)
    println("done.")

end

# ----------------- SIMULATION OPTIONS -----------------

out_path = joinpath("output","run_"*Dates.format(Dates.now(), "yyyy_mm_dd__HH_MM_SS"))
mkdir(out_path)

# FROG data file
FROG_dir = joinpath("input", "tangerine_FROG_50kHz")

file_FROG_Et = "Ek.dat"
file_FROG_Eω = "Speck.dat"

path_FROG_Eω = joinpath(FROG_dir, file_FROG_Eω)

# mirror data files
path_m1 = joinpath("input", "mirrordata_0fs2.txt")
path_m2 = joinpath("input", "mirrordata_-30fs2.txt")

#comments = "Simulation with reflectivity from mirror data file but GDD set constant (0 and -30.0 fs2)." #constant 100% reflectivity but GDD taken from mirror data file."
comments = " "

# ----------------- SETTING RESULT PROCESSING -----------------
use_FROG_spectrum = false
show = false                    # if true: show plots
show_title = true
save = true
use_pdf = false
norm = true
disable_latex = true # if true: disable latex rendering of plots (saves time but might result in some labels or title being displayed incorrectly)


@eval import PyPlot: pygui, plt, PyDict, matplotlib 
    close("all")
    pygui(true)

    # set plot formatting 
    rcParams = PyDict(matplotlib."rcParams") # get rcParams 
    if disable_latex==false  rcParams["text.usetex"] = true end # enable LaTeX renadering
    rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
    rcParams["font.family"] = "STIXGeneral" # use LateX font for text
    rcParams["font.size"] = 16 # set standard font size 
    fig_dim = 2* [3.14961, 2.3622075] # for 8cm width ; double for 16cm width  


    #plt.rc("text", usetex=false)


# ----------------- SETTING PARAMETERS -----------------

# pulse parameters
λ0 = 1030e-9        # central wavelength [m]
τ = 150e-15         # pulse duration [s]; ignored if FROG spectrum is used
E_pulse = 250e-6    # pulse energy [J]
w0 = 150e-6         # beam waist at focus [m]
#M2 = 1.16           # beam quality     # not yet implemented; how would you do that?

# MPC parameters
L = 380.1934e-3          # MPC cell length [m]
R = 200e-3          # radius of curvature [m]
propz = -L/2         # propagation distance from the waist [m] (distance between mirror position and focal position)
z_vals = L .* [0, 1/2, 1]     # points along the cell at which to investigate beam evolution [m]

N_reflections = 1
Npass = 30

pres = 1.5               # gas pressure [bar]
gas = :Kr           # gas type

ion = true         # if true: enable ionisation response, if false: disable ionisation 
ion_model="PPT"    # set to "ADK" or "PPT" (has no effect if ion==false); "ADK" is less accurate at low intensities but faster; "PPT" may crash at very high intensities


# simulation grid parameters
λlims = (700e-9, 1300e-9) #(600e-9, 1500e-9) # wavelength range of interest [m]
trange = 20*τ #300.0e-15 #0.05e-12    # total extent of time window required [s], default was 50fs (NOTE: if this is too short, the range is extended automatically; this value was taken from an older version of the code; no justification!)
Nz = 201                            # number of points along z at which the spectrum is saved


# calculate theoretical beam size in the MPC (as in https://doi.org/10.1364/OPTICA.449225 eq. 3 and 4)
    C = L/R 

    # beam waist in the focus
    w02 = R*λ0 * sqrt(C*(2-C)) /2pi
    w0 = sqrt(w02)                                                # ! this value is used to define the gaussian beam in the beginning of the simulation
    @printf("Input focal spot radius w0 = %.2f µm\n", w0*1e6)

    # beam waist on the mirror
    wm2 = R*λ0 * sqrt(C/(2-C)) /pi
    wm = sqrt(wm2) 
    @printf("Calculated spot radius at the mirror location wm = %.2f µm\n", wm*1e6)


# initialize arrays to track beam size at L/2 and L (in the center of the cell and at the mirror position)
    w_0 = []
    w_m = []


# Hankel transformation
R_hankel = 3.0*wm           # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  
N_hankel = 256 #512 #1024             # sample size for Hankel tansform grid           


# ----------------- SET SIMULATION GRID ----------------------------

q = Hankel.QDHT(R_hankel, N_hankel, dim=2)                  # set up discrete Hankel transform matrix 
grid = Grid.RealGrid(L, λ0, λlims, trange)               # set up time & space grid for gradient approximation 
                   
energyfun, energyfun_ω = Fields.energyfuncs(grid, q)    # "energyfun" gives total energy in a field E(t); energyfun_ω is needed in case field is defined via FROG spectrum  


# ----------------- SET NONLINEAR EFFECTS ----------------------------

ionpot = PhysData.ionisation_potential(gas)                 # set gas ionisation potential   

if ion_model=="ADK"
    ionrate = Ionisation.ionrate_fun!_ADK(gas)                  # set gas ionisation rate (ADK)
elseif ion_model=="PPT"
    ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)        # set gas ionisation rate (PPT)
end    

n_gas = PhysData.ref_index_fun(gas, pres)                       # refractive index of the gas for given pressure (no pressure profile, pressure is constant)

linop = LinearOps.make_const_linop(grid, q, n_gas)              # generate linear operator for pulse-propagation equation
normfun = NonlinearRHS.const_norm_radial(grid, q, n_gas)        # generate normalisation function for radial symmetry


# * * * SET KERR EFFECT AND/OR PLASMA FORMATION 
if ion == true                             # nonlinear response function with ionisation and Kerr effect
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
            Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),)  
elseif ion == false                        # nonlinear response function without ionisation, just Kerr effect
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)    
end 

# ---------------- SET GAS DENSITY ----------------
densityfun = let dens0=PhysData.density(gas, pres)      # pressure is constant in MPC
    z -> dens0 
end

# ---------------- SET UP SIMULATION -------------------
if use_FROG_spectrum==true
    inputs = () # leave input field empty as we will be defined later by the FROG spectrum

else
    # input field
    # to get a converging beam, we define a Gaussian spot with size w0
    # and then *back-propagate* by L/2 i.e. by half of our propagation window
    # this means the nominal focus is halfway through our propagation in this example
    inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ, energy=E_pulse, w0=w0, propz=propz)
end

Eωk, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)  # set up propagation; Eωk is overwritten when FROG spectrum is used


# ---------------- DEFINE MIRRORS ----------------

# MIRROR 1
    data = readdlm(path_m1, skipstart=1)

    λ_m1 = data[:,1].*1e-9     # read in wavelengths [m]
    R_m1 = data[:,2].*1e-2     # read in reflectivity [fractions]
    GDD_m1 = data[:,3].*1e-30  # read in GDD [s]


# MIRROR 2
    data = readdlm(path_m2, skipstart=1)

    λ_m2 = data[:,1].*1e-9     # read in wavelengths [m]
    R_m2 = data[:,2].*1e-2     # read in reflectivity [fractions]
    GDD_m2 = data[:,3].*1e-30  # read in GDD [s]



# ----------------- DEFINE INPUT VIA COMPLEX FROG SPECTRUM -----------------
if use_FROG_spectrum==true

    println("using FROG spectrum as input")

    # READ IN AND PREPARE FROG DATA 
        data = readdlm(path_FROG_Eω)

        λdat = data[:,1]     # read in wavelengths [nm]
        indices = findall(x -> x >= 0, λdat)     # cut off negative wavelengths
        λdat = λdat[indices]

        ωdat = PhysData.wlfreq.(λdat*1e-9)  # calculate angular frequencies 

        E_real = data[indices,4]   # read in real part of complex field
        E_imag = data[indices,5]   # read in imaginary part of complex field

        Eωdat = E_real .+ 1im.*E_imag   # combine real and imaginary part
        println("size(Eω_FROG): ", size(Eωdat))


    # CREATE DATAFIELD (see https://github.com/LupoLab/Luna.jl/issues/357)
        df = Fields.DataField(ωdat, Eωdat; energy=1) # energy is irrelevant here as we will rescale again
        Eω = df(grid, nothing) # second argument is unused for DataField (need the planned FT for other types)

        Eωr = Eω .* sqrt.(Maths.gauss.(q.r, w0/2))' # combine with spatial distribution; beam waist taken from above calculation

        Eωk = q * Eωr # transform from ω-r space to ω-k space
        Eωk .*= sqrt(E_pulse)/sqrt(energyfun_ω(Eωk)) # rescale to actual energy

        Fields.prop!(Eωk, propz, grid, q)           # propagate linearly back from the waist to the position of the mirror where simulation starts
end


println(size(Eωk))
Eωk_init = copy(Eωk)


# transform to E(ω, r)
Eωr = q \ Eωk # inverse QDHT: Eωk is E(ω, k) so Eωr is E(ω, r)


# apply reflectivity and GDD of the MPC mirror  
for pass in 1:Npass
    if isodd(pass)
        reflecting_mirror1!(Eωr)
    else
        reflecting_mirror2!(Eωr) 
    end
end


# transform back to E(ω, k)
Eωk .= q*Eωr        # the '.' ensures that the input is updated and not assigned to a new variable

# time domain after mirror
Ẽω_init = zeros(ComplexF64, size(Eωk_init, 1))      # set up new array for Ẽω = Ẽ(ω, z); COMPLEX electric field amplitude in FREQUENCY domain INTEGRATED along r (from 0 to infinity to my knowledge)
Ẽω = zeros(ComplexF64, size(Eωk, 1))      # set up new array for Ẽω = Ẽ(ω, z); COMPLEX electric field amplitude in FREQUENCY domain INTEGRATED along r (from 0 to infinity to my knowledge)


for i = 1:size(Eωk_init, 1)
    Ẽω_init[i] = Hankel.integrateK(Eωk_init[i,:], q)
end
Iω_init = abs2.(Ẽω_init)
ϕω_init = -angle.(Ẽω_init)

for i = 1:size(Eωk, 1)
    Ẽω[i] = Hankel.integrateK(Eωk[i,:], q)
end
Iω = abs2.(Ẽω)
ϕω = -angle.(Ẽω)

Et_init = FFTW.irfft(Ẽω_init, length(grid.t),1)
It_init = abs2.(Maths.hilbert(Et_init))
t_peak_init = grid.t[argmax(It_init)]                             # time for maximum intensity [s]
ceil(t_peak_init, digits=1)

Et = FFTW.irfft(Ẽω, length(grid.t),1)
It = abs2.(Maths.hilbert(Et))
t_peak = grid.t[argmax(It)]                       # time for maximum intensity [s]

λ = PhysData.wlfreq.(grid.ω)


#+++++ PLOT 1:  Frequency domain: Spectrum before and after a number of passes
fig, ax1 = plt.subplots(figsize=fig_dim)
if show_title ax1.set_title("Spectrum before and after $(Npass) pass(es)") end
ax1.plot(λ[2:end]*1e9, norm ? Iω_init[2:end]/maximum(Iω_init[2:end]) : Iω_init[2:end], label="before reflection")#, color="red")
ax1.plot(λ[2:end]*1e9, norm ? Iω[2:end]/maximum(Iω[2:end]) : Iω[2:end], label="after reflection", linestyle="--")#, color="red")
ax2 = ax1.twinx()
ax2.plot(λ[2:end]*1e9, unwrap(ϕω_init[2:end]), label="phase")#, color="orange")
ax2.plot(λ[2:end]*1e9, unwrap(ϕω[2:end]), label="phase", linestyle="--")#, color="orange")

ax1.set_xlim(960, 1100)
#ax1.set_xlim(λlims[1]*1e9, λlims[2]*1e9)
#ax2.set_ylim(-50, 50)
ax1.set_xlabel("λ [nm]")
ax1.set_ylabel(norm ? "I (norm.)" : "I (arb. units)")
ax2.set_ylabel("phase [rad]")
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
plt.tight_layout()

#+++++ PLOT 1B:  Frequency domain: Spectrum before and after reflection(s)
fig, ax1 = plt.subplots(figsize=fig_dim)
if show_title ax1.set_title("Spectrum before and after $(Npass) pass(es)") end
ax1.plot(grid.ω.*1e-15, norm ? Iω_init/maximum(Iω_init) : Iω_init, label="before reflection")#, color="red")
ax1.plot(grid.ω.*1e-15, norm ? Iω/maximum(Iω) : Iω, label="after reflection", linestyle="--")#, color="red")
ax2 = ax1.twinx()
ax2.plot(grid.ω.*1e-15, unwrap(ϕω_init), label="phase")#, color="orange")
ax2.plot(grid.ω.*1e-15, unwrap(ϕω), label="phase", linestyle="--")#, color="orange")

ax1.set_xlim(1.7, 1.95)
#ax1.set_xlim(λlims[1]*1e9, λlims[2]*1e9)
#ax2.set_ylim(-50, 50)
ax1.set_xlabel("ω [PHz]")
ax1.set_ylabel(norm ? "I (norm.)" : "I (arb. units)")
ax2.set_ylabel("phase [rad]")
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
plt.tight_layout()

#+++++ PLOT 2:  Time domain: Pulse before and after refection(s)
plt.figure(figsize=fig_dim) 
if show_title plt.title("Pulse before and after $(Npass) pass(es)") end
plt.xlabel("t (fs)")
plt.xlim(minimum(grid.t)*1e15, maximum(grid.t)*1e15)
plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")

# create plot label 1
    formatted_number = @sprintf("%.2f", t_peak_init.*1e15)
    label_text = "before reflection \nt_peak = $(formatted_number) fs"
plt.plot(grid.t*1e15, norm ? Maths.normbymax(It_init[:,1]) : It_init[:,1], label=label_text)#, label=L"\tau_{FWHM}="*string(round(τ_input, digits=1) )*"fs")

# create plot label 2
    formatted_number = @sprintf("%.2f", t_peak.*1e15)
    label_text = "after reflection \nt_peak = $(formatted_number) fs"
plt.plot(grid.t*1e15, norm ? Maths.normbymax(It[:,1]) : It[:,1], label=label_text, linestyle="--")#, label=L"\tau_{FWHM}="*string(round(τ_input, digits=1) )*"fs")

plt.legend(loc="upper right")
plt.tight_layout()

plt.show()

# if save==true
#     if use_pdf == true 
#         plt.savefig(joinpath(out_path,"time_domain_pass_$(pass)_after_reflection.pdf"))
#     else     
#         plt.savefig(joinpath(out_path,"time_domain_pass_$(pass)_after_reflection.png"),dpi=1000)
#     end 
# end

# if show == false
#     plt.close()
# end
