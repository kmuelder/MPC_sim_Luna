using Luna
using LaTeXStrings

@eval import PyPlot: pygui, plt, PyDict, matplotlib 
    close("all")
    pygui(true)

    # # set plot formatting 
    # rcParams = PyDict(matplotlib."rcParams") # get rcParams 
    # if disable_latex==false  rcParams["text.usetex"] = true end # enable LaTeX renadering
    # rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
    # rcParams["font.family"] = "STIXGeneral" # use LateX font for text
    # rcParams["font.size"] = 16 # set standard font size 
    # fig_dim = 2* [3.14961, 2.3622075] # for 8cm width ; double for 16cm width
"""
To run a simple simulation of ultrafast pulse propagation in a gas-filled hollow capillary fibre, 
you can use prop_capillary. As an example, take a 3-metre length of HCF with 125 μm core radius, 
filled with 1 bar of helium gas, and driving pulses centred at 800 nm wavelength with 120 μJ of energy and 10 fs duration. 
We consider a frequency grid which spans from 120 nm to 4 μm and a time window of 1 ps.

output = prop_capillary(125e-6, 3, :He, 1; λ0=800e-9, energy=120e-6, τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12)

Plotting.prop_2D(output)

You can also display the power spectrum at the input and output (and anywhere in between):
Plotting.spec_1D(output, [0, 1.5, 3]; log10=true)
"""

R = 300e-3          # radius of curvature in mm
L = 290e-3          # MPC cell length in mm
Nrt = 35            # number of round trips
p = 1               # gas pressure in bar
gas = :Xe           # gas type
λ0 = 1030e-9        # central wavelength in nm
λlims = (950e-9, 1100e-9) # wavelength range of interest (?) in nm
n2 = 5.2e-23        # nonlinear index in m^2/W
τ = 300e-15         # pulse duration in fs
E_pulse = 150e-6    # pulse energy in microJoule
A_eff = 0.19e-6     # effective area of equivalent waveguide in mm^2
P_ratio = 0.18      # ratio of peak power over critical power (=critical power for self-focusing)

N_grid = 256        # sampling size of time axis (?)
dt = 10e-15         # time step in fs
trange = N_grid*dt  # time window in fs

L_tot = 2*Nrt*L       # equivalent waveguide length
println("total propagation length = ",L_tot)

M = 1
A_eff_cal = M^2 * λ0*L/(2*atan(sqrt(L/(2*R-L))))
println("A_eff_cal = ", A_eff_cal, "[m^2]")

#r = sqrt(2*A_eff_cal/pi)          # fiber radius in micrometer
r = 360.2e-6 #400e-6 #360.2e-6

println("fiber radius r = ", r*1e6, " [microns]")


### CALCULATING FIBER PROPAGATION
output = prop_capillary(r, L_tot, gas, p; λ0=λ0, energy=E_pulse, τfwhm=τ, λlims=λlims, trange=trange)

### PLOTTING SPECTRAL AND TEMPORAL EVOLUTION WHILE PROPAGATION
#Plotting.prop_2D(output, trange=(-2560e-15, 2560e-15), λrange=λlims)

### PLOTTING SPECTRUM AT DIFFERENT POINTS OF PROPAGATION
spectrum_1D = Plotting.spec_1D(output, [L_tot]; λrange=λlims, log10=false)



"""
In this section I'm trying to find the core radius for which the effective has the same value as given in the paper.
I am using a function from luna that calculates the effective area. This function needs an AbstractMode type as input.
So far I only found the so called MarcatiliModes as Abstract Mode types in luna.
I am creating such a mode "for a capillary with radius `a` which is filled with `gas` to pressure `P`".
The effective area is calculated for a range of core radius values.
"""

r = 350e-6: 1e-6: 370e-6
gas = :Xe
P = 1
area = zeros(length(r))

for (i, a) in enumerate(r)
    mode = Capillary.MarcatiliMode(a, gas, P)
    A_eff = Modes.Aeff(mode)
    println("Aeff(r = $(a) [m]) = ", A_eff)
    area[i] = A_eff
end

plt.figure()
plt.title("calculated effective mode areas for different core radii")
plt.scatter(r.*1e6, area)
plt.axhline(1.9428e-7, color="black", linestyle="dashed", label="Aeff from paper")
plt.xlabel("core radius r [microns]")
plt.ylabel("effective mode area [m^2]")
plt.legend()
plt.show()
