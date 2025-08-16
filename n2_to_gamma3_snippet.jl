n2_data = 5.5e-23 # data from Daher et al for Xenon
pressure_data = 1 # 1 bar
density_data = PhysData.density(gas, pressure_data)
χ3 = 4/3 * PhysData.ε_0*PhysData.c * n2_data # note: assuming n0 ≈ 1
γ3 = χ3/density_data

responses = (Nonlinear.Kerr_field(γ3),
             Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))