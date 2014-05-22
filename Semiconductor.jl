
module Semiconductor

# Teperature dependence of energy gap 
# varshni(e0, a, b, t)
export varshni
function varshni(e0::Float64, a::Float64, b::Float64, t::Float64) 
    return e0 - a * t * t / (t + b)
end

# Vegard's rule to approximate alloys' properties with bowing term
export vegard
function vegard(x::Float64, param1::Float64, param2::Float64, bow::Float64) 
	if !(0.0 <= x <= 1.0)
		warn("vegard: x is expected to be 0<=x<=1")
	end
    return param1 * x + param2 * (1 - x) - bow * x * (1 - x)
end

end