using Luna
using LaTeXStrings

@eval import PyPlot: pygui, plt, PyDict, matplotlib 
    close("all")
    pygui(true)

"""
FROM luna's READ ME on the prop_capillary function: 
    To run a simple simulation of ultrafast pulse propagation in a gas-filled hollow capillary fibre, 
    you can use prop_capillary. As an example, take a 3-metre length of HCF with 125 μm core radius, 
    filled with 1 bar of helium gas, and driving pulses centred at 800 nm wavelength with 120 μJ of energy and 10 fs duration. 
    We consider a frequency grid which spans from 120 nm to 4 μm and a time window of 1 ps.

    output = prop_capillary(125e-6, 3, :He, 1; λ0=800e-9, energy=120e-6, τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12)

    Plotting.prop_2D(output)

    You can also display the power spectrum at the input and output (and anywhere in between):
    Plotting.spec_1D(output, [0, 1.5, 3]; log10=true)
"""


figure_1 = true
figure_2 = true
figure_5 = true
esmerando_benchmarking = true



### SELLMEIER RELATION CITED IN THE PAPER TO CALCULATE βs ###
    """
    Sellmeier relation that calculates the refractive index n(λ) for Xenon taken from:
        A. Bideau-Mehu, Y. Guern, R. Abjean, and A. Johannin-Gilles, 
        “Measurement of refractive indices of neon, argon, krypton, and xenon in the 253.7-140.4 nm wavelength range. 
        Dispersion relations and estimated oscillator strengths of the resonance lines,” 
        J. Quant. Spectrosc. Radiat. Transfer 25, 395–402 (1981).
    which is used in MPC/fiber paper. 
    """
    function sellmeier_Xe(λ)
        λ *= 1e6 # converting to micrometer
        n = @. 1.2055e-2*(0.26783/(46.301-(1/λ^2)) + 0.29481/(50.578-(1/λ^2)) + 5.0333/(112.74-(1/λ^2))) + 1 # not sure if its right to just ad +1 here, the formula without that calculates (n-1)
        return n 
    end

    λ = 120e-9:1e-9:1100e-9 # wavelength in [m]

    plt.figure()
    plt.title("Dispersion formula for Xenon")
    plt.plot(λ.*1e9, sellmeier_Xe(λ).*1e6 .-1)
    plt.xlabel("λ [nm]")
    plt.ylabel("(n-1)")
    plt.xlim(120, 260)
    #plt.ylim(0, 6000)
    plt.show()



### SIMULATION FOR FIG.1 FROM 1D MPC/FIBER PAPER: https://doi.org/10.1364/JOSAB.386049 ###
if figure_1 == true
    # SETTING PARAMETERS

        # pulse parameters
            λ0 = 1030e-9        # central wavelength in nm
            τ = 300e-15         # pulse duration in fs
            E_pulse = 150e-6    # pulse energy in microJoule
            n2 = 5.2e-23        # nonlinear index in [m^2/W]
            M2 = 1.0               # beam quality

        # MPC parameters
            L = 290e-3          # MPC cell length in mm
            R = 300e-3          # radius of curvature in mm
            Nrt = 35            # number of round trips

            pres = 1               # gas pressure in bar
            gas = :Xe           # gas type

        
        # simulation grid parameters
            λlims = (750e-9, 1300e-9) #(950e-9, 1100e-9) # wavelength range of interest (?) in nm
            
            N_grid = 1024        # sampling size of time axis (?)
            dt = 5e-15         # time step in fs
            trange = N_grid*dt  # time window in fs
            


    # CALCULATE INPUT FOR LUNA FUNCTION

        flength = 2*Nrt*L     # equivalent waveguide length

        A_eff = M2*λ0*L / (2*atan(sqrt(L/(2R-L))))     # effective area of equivalent waveguide in [m^2]
        println("A_eff = ",A_eff)

        k0 = 2π/λ0          # wavenumber of central wavelength [rad/m]
        γ = k0*n2/A_eff     # nonlinear coefficient in [rad/(m*W)]
        println("γ = ",γ)


    ### SIMULATION DONE WITH β FROM LUNA ###

            # Calculate β for given gas & pressure
                orders = 0:1:7
                βs = zeros(length(orders)) 

                for (i, order) in enumerate(orders)
                    βs[i] = PhysData.dispersion(order, gas, λ0, pres) 
                end
                println("βs luna: ",βs)



        # RUN SIMULATION

            output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)


        # PLOT RESULTS

            Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(950e-9, 1100e-9), trange=(-1e-12, 1e-12))

            #zslice = [0, flength/2, flength]
            zslice = [flength]
            Plotting.spec_1D(output, zslice; λrange=(950e-9, 1100e-9), log10=false)
            plt.title("Using lunas βs values")
            plt.tight_layout()


            Plotting.time_1D(output, zslice; trange=(-500e-15, 500e-15))#, bandpass=(180e-9, 220e-9))

            sellmeier_coeff = PhysData.sellmeier_gas(:Xe)
            println(sellmeier_coeff)



    ### SIMULATION DONE WITH β THEY USED IN THE MPC/FIBER PAPER ###

        # Calculate β for given gas & pressure
            orders = 0:1:7
            βs = zeros(length(orders)) 

            for (i, order) in enumerate(orders)
                βs[i] = PhysData.dispersion_func(order, sellmeier_Xe).(λ0)  
            end
            println("βs manually: ", βs)


        # RUN SIMULATION

        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)


        # PLOT RESULTS

            Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(950e-9, 1100e-9), trange=(-1e-12, 1e-12))

            #zslice = [0, flength/2, flength]
            zslice = [flength]
            Plotting.spec_1D(output, zslice; λrange=(950e-9, 1100e-9), log10=false)
            plt.title("Using papers βs values")
            plt.tight_layout()


            Plotting.time_1D(output, zslice; trange=(-500e-15, 500e-15))#, bandpass=(180e-9, 220e-9))

            sellmeier_coeff = PhysData.sellmeier_gas(:Xe)
            println(sellmeier_coeff)

end


### SIMULATION FOR FIG.2 FROM 1D MPC/FIBER PAPER: https://doi.org/10.1364/JOSAB.386049 ###
if figure_2 == true
    # SETTING PARAMETERS

        # pulse parameters
        λ0 = 1030e-9        # central wavelength in nm
        τ = 350e-15         # pulse duration in fs
        E_pulse = 200e-6    # pulse energy in microJoule
        
        M2 = 1.0               # beam quality

    # MPC parameters
        L = 551e-3          # MPC cell length in mm
        R = 300e-3          # radius of curvature in mm
        Nrt = 5            # number of round trips

        pres = 3.0               # gas pressure in bar
        gas = :Xe           # gas type
        use_paper_βs = true

    
    # simulation grid parameters
        λlims = (750e-9, 1300e-9) #(950e-9, 1100e-9) # wavelength range of interest (?) in nm
        
        N_grid = 1024        # sampling size of time axis (?)
        dt = 5e-15         # time step in fs
        trange = N_grid*dt  # time window in fs
        


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


    if use_paper_βs == false
    # Calculate β for given gas & pressure from lunas data
        orders = 0:1:7
        βs = zeros(length(orders)) 

        for (i, order) in enumerate(orders)
            βs[i] = PhysData.dispersion(order, gas, λ0, pres) 
        end
        println("βs luna: ",βs)

    elseif use_paper_βs == true
    # Calculate β from Sellmeier relation cited in paper
        orders = 0:1:7
        βs = zeros(length(orders)) 

        for (i, order) in enumerate(orders)
            βs[i] = PhysData.dispersion_func(order, sellmeier_Xe).(λ0)  
        end
        println("βs manually: ", βs)
    else
        println("use_paper_βs has to be either true or false")
    end



    # RUN SIMULATION

        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)


    # PLOT RESULTS

        Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(950e-9, 1100e-9), trange=(-1e-12, 1e-12))

        #zslice = [0, flength/2, flength]
        zslice = [flength]
        Plotting.spec_1D(output, zslice; λrange=(950e-9, 1100e-9), log10=false)
        plt.title(use_paper_βs ? "Using papers βs" : "Using lunas βs values")
        plt.tight_layout()


        Plotting.time_1D(output, zslice; trange=(-500e-15, 500e-15))#, bandpass=(180e-9, 220e-9))

        sellmeier_coeff = PhysData.sellmeier_gas(:Xe)
        println(sellmeier_coeff)

end

### SIMULATION FOR FIG.5 FROM 1D MPC/FIBER PAPER: https://doi.org/10.1364/JOSAB.386049 ###
if figure_5 == true

    # SETTING PARAMETERS

        # pulse parameters
        λ0 = 1030e-9        # central wavelength in nm
        τ = 338e-15         # pulse duration in fs
        E_pulse = 33e-6    # pulse energy in microJoule
        
        M2 = 1.2*1.2               # beam quality

        # MPC parameters
        L = 551e-3          # MPC cell length in mm
        R = 300e-3          # radius of curvature in mm
        Nrt = 21            # number of round trips

        GDD = -25.0e-30       # GDD per pass due to chirped mirrors in fs^2

        pres = 3.8               # gas pressure in bar
        gas = :Xe           # gas type
        use_paper_βs = false

    
        # simulation grid parameters
        λlims = (750e-9, 1300e-9) #(950e-9, 1100e-9) # wavelength range of interest (?) in nm
        
        N_grid = 1024        # sampling size of time axis (?)
        dt = 5e-15         # time step in fs
        trange = N_grid*dt  # time window in fs
        


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


    if use_paper_βs == false
    # Calculate β for given gas & pressure from lunas data
        orders = 0:1:7
        βs = zeros(length(orders)) 

        for (i, order) in enumerate(orders)
            βs[i] = PhysData.dispersion(order, gas, λ0, pres) 
        end
        println("βs luna: ",βs)

        # add mirror GDD
        βs[3] += GDD/L
        println("βs luna chirped: ",βs)

    elseif use_paper_βs == true
    # Calculate β from Sellmeier relation cited in paper
        orders = 0:1:7
        βs = zeros(length(orders)) 

        for (i, order) in enumerate(orders)
            βs[i] = PhysData.dispersion_func(order, sellmeier_Xe).(λ0)  
        end
        println("βs paper: ", βs)

        # add mirror GDD
        βs[3] += GDD/L
        println("βs paper chirped: ",βs)
    else
        println("use_paper_βs has to be either true or false")
    end



    # RUN SIMULATION

        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)


    # PLOT RESULTS

        Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(930e-9, 1130e-9), trange=(-1e-12, 1e-12))

        #zslice = [0, flength/2, flength]
        zslice = [flength]
        Plotting.spec_1D(output, zslice; λrange=(930e-9, 1130e-9), log10=false)
        plt.title(use_paper_βs ? "Using papers βs" : "Using lunas βs values")
        plt.tight_layout()


        Plotting.time_1D(output, zslice; trange=(-500e-15, 500e-15))#, bandpass=(180e-9, 220e-9))

        sellmeier_coeff = PhysData.sellmeier_gas(:Xe)
        println(sellmeier_coeff)

end

# ### SIMULATION FOR https://doi.org/10.1364/OL.43.005877 AS ESMERANDO COMPARED HIS RESULTS TO THAT ###
if esmerando_benchmarking == true

    # SETTING PARAMETERS

        # pulse parameters
        λ0 = 1030e-9        # central wavelength in nm
        τ = 1.3e-12         # pulse duration 
        E_pulse = 18.6e-3    # pulse energy 

        M2 = 1.17*1.18 #1.10*1.09 #1.17*1.18 #             # beam quality

    # MPC parameters
        L = 1.99*1.5          # MPC cell length in m
        R = 1.5          # radius of curvature in m
        Nrt = 45/2            # number of round trips (in the paper they give the number of passes)

        pres = 600e-3               # gas pressure in bar
        gas = :Ar           # gas type

        


    # simulation grid parameters
        λlims = (700e-9, 1500e-9)#(550e-9, 1500e-9) #(950e-9, 1100e-9) # wavelength range of interest (?) in nm
        
        N_grid = 2048        # sampling size of time axis (?)
        dt = 1.0e-15         # time step in ps
        trange = 7*τ #10.0e-12 #N_grid*dt  # time window in fs
        


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


    ### SIMULATION DONE WITH β FROM LUNA ###

        # Calculate β for given gas & pressure
            orders = 0:1:7
            βs = zeros(length(orders)) 

            for (i, order) in enumerate(orders)
                βs[i] = PhysData.dispersion(order, gas, λ0, pres) 
            end
            println("βs luna: ",βs)

            

    # RUN SIMULATION

        output = prop_gnlse(γ, flength, βs; λ0, λlims, trange, τfwhm=τ, energy=E_pulse, pulseshape=:gauss)


    # PLOT RESULTS

        Plotting.prop_2D(output, :λ, dBmin=-40.0,  λrange=(980e-9, 1070e-9), trange=(-5e-12, 5e-12))

        #zslice = [0, flength/2, flength]
        zslice = [flength]
        Plotting.spec_1D(output, zslice; λrange=(980e-9, 1070e-9), log10=false)
        plt.title("M2(input) = $(ceil(M2 * 100) / 100)")
        plt.tight_layout()


        Plotting.time_1D(output, zslice; trange=(-5e-12, 5e-12))#, bandpass=(180e-9, 220e-9))
end