module constants

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constants

export CONST_PLANCK, CONST_PLANCK_REDUCED
const  CONST_PLANCK 		= 6.6260695729e-27  				# erg * sec
const  CONST_PLANCK_REDUCED = CONST_PLANCK / 2pi  				# erg * sec

export ELECTRON_CHARGE, ELECTRON_CHARGE_SI, ELECTRON_MASS
const  ELECTRON_CHARGE 		= 4.803529695e-10   				# CGS
const  ELECTRON_CHARGE_SI 	= 1.60217656535e-19 				# Coloumbs
const  ELECTRON_MASS 		= 9.10938291e-28    				# grams

export CONST_BOLTZMANN
const  CONST_BOLTZMANN 		= 1.3806488e-16     				# erg / K

export SPEED_LIGHT
const  SPEED_LIGHT 			= 29979245800    					# cm per sec

export EV2NM, NM2EV
const  EV2NM = CONST_PLANCK * SPEED_LIGHT / ELECTRON_CHARGE_SI  # eV * nm
const  NM2EV = EV2NM 


end