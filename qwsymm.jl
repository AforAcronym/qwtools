#!/usr/bin/env julia

# tic()

using constants
using utils



function qwsymm( V::Real, width::Real, mass_w::Real, mass_b::Real )
# qwsymm       gives quantum levels of a single quantum well (QW) with 
#              symmetric barriers
#
# qwsymm( V, width, mass_w, mass_b )
#
# V         -- barriers width (depth of the QW) in ergs or electronvolts
# width     -- width of the QW in cm
# mass_w    -- a particle mass in the QW either in units of free electron
#              mass or in grams
# mass_b    -- a particle mass in the barriers either in units of free 
#              electron mass or in grams
#
	# converting units
	if mass_w > ELECTRON_MASS * 15
		mass_w = mass_w * ELECTRON_MASS 
	end

	if mass_b > ELECTRON_MASS * 15
		mass_b = mass_b * ELECTRON_MASS 
	end

	if V > 0.001
		V = convert_ev2erg(V) 
	end

	# Calculations
	# TODO
	wavenum_w_xy = 0 

	# Wavenumber in the well
	wavenum_w_z(E::Real) 	= sqrt( 2 * mass_w *    E    / CONST_PLANCK_REDUCED^2 - wavenum_w_xy^2 )

	# Wavenumber in barriers
	wavenum_b_z(E::Real) 	= sqrt( 2 * mass_b * (V - E) / CONST_PLANCK_REDUCED^2 + wavenum_w_xy^2 )
	etha(E::Real) 		 	= mass_w / mass_b * (wavenum_b_z(E) / wavenum_w_z(E))

	# Curves for finding roots
	solution_even(E::Real) 	= tan(wavenum_w_z(E) * width / 2) - etha(E) 
	solution_odd(E::Real)  	= cot(wavenum_w_z(E) * width / 2) + etha(E) 

	tol = eps(Float64); # tolerance in erg, about 6e-05 ev
	energyeigen = [];
	err = [];

	# Step is 0.1 meV, desired initial resolution
	energystep = convert_ev2erg(1e-4)
	energy = 0
	energynext = energystep

	while energy < V
		
		if (energynext > V)
			# Avoid complex solutions that lay above the barriers
			energynext = V
		end

		if solution_even(energy) * solution_even(energynext) < 0  &&  solution_even(energy) < 0
			energy_level, level_err = bisection(solution_even, energy, energynext, tol);
			energyeigen = [energyeigen,  energy_level ];
			err = [err, level_err ];
		end

		if solution_odd(energy) * solution_odd(energynext)  < 0  &&  solution_odd(energy) < 0
			energy_level, level_err = bisection(solution_odd, energy, energynext, tol);
			energyeigen = [energyeigen,  energy_level ];
			err = [err,  level_err ];
		end

		energy = energynext
		energynext += energystep
	end

	return (energyeigen, err)
end






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test qwsymm function 
# qwsymm(0.4, 3 * 1e-7, 0.067, 0.067) must result in 0.178 eV
# and print 0.17826419124082454 ± 4.7267354010182736e-17 eV

# Get value and error arrays, convert ergs to eV, print the result
function qwsymm_ev_print( V::Real, width::Real, mass_w::Real, mass_b::Real )
	energy_erg, err_erg = qwsymm( V, width, mass_w, mass_b )
	energy_ev 	= convert_erg2ev(energy_erg)
	err_ev 		= convert_erg2ev(err_erg)
	println( join( [string(energy_ev[i]) * " ± " * string(err_ev[i]) * " eV" for i = 1:length(energy_ev)] ,"\n") )

end


# toc()
