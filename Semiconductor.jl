
module Semiconductor


# ---------------------------------------------------
# Teperature dependence of energy gap 
# varshni(e0, a, b, t)
export varshni
function varshni(e0::Float64, a::Float64, b::Float64, t::Float64) 
    return e0 - a * t * t / (t + b)
end


# ---------------------------------------------------
# Vegard's rule to approximate alloys' properties with bowing term
export vegard
function vegard(x::Float64, param1::Float64, param2::Float64, bow::Float64) 
	if !(0.0 <= x <= 1.0)
		warn("vegard: x is expected to be 0<=x<=1")
	end
    return param1 * x + param2 * (1 - x) - bow * x * (1 - x)
end

# General crystal lattice
abstract Lattice



# ---------------------------------------------------
type Cubic <: Lattice 
	a::Float64  
end

function cellvolume(l::Cubic)
	return l.a * l.a * l.a
end



# ---------------------------------------------------
type Wurzite <: Lattice 
	a::Float64  
	c::Float64 
end

function cellvolume(l::Wurzite)
	return l.c * l.a * l.a * 1.5 * sqrt(3)
end


# ---------------------------------------------------
type Semiconductor
	energy_gap_0::Float64
	lattice::Lattice
	polarization::Vector{Float64,1}
	piezotensor::Tensor #?
end

end