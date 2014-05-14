module utils

using constants

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export convert_ev2nm, convert_nm2ev, convert_ev2erg, convert_erg2ev

# Convert from electrnvolts to nanometers
convert_ev2nm(energy::Real) = NM2EV / energy
convert_ev2nm{T<:Real}(energy::Array{T}) = NM2EV / energy

# Convert from nanometers to electrnvolts
convert_nm2ev(wavelength::Real) = NM2EV / wavelength
convert_nm2ev{T<:Real}(wavelength::Array{T}) = NM2EV / wavelength

# Convert from electrnvolts to ergs (CGS)
convert_ev2erg(ev::Real) = ELECTRON_CHARGE_SI * 1e7 * ev
convert_ev2erg{T<:Real}(ev::Array{T}) = ELECTRON_CHARGE_SI * 1e7 * ev

# Convert from ergs (CGS) to electrnvolts
convert_erg2ev(erg::Real) = erg / (ELECTRON_CHARGE_SI * 1e7)
convert_erg2ev{T<:Real}(erg::Array{T}) = erg / (ELECTRON_CHARGE_SI * 1e7)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export bisection

# My implementation of the bisection method
function bisection(func::Function, a::Real, b::Real, tol::Real)
	if (b < a)
		a, b = b, a
	end

	# Check for positive sign
	tol = tol < 0 ? -tol : tol  
	tol = tol < eps(Float64) ? eps(Float64) : tol

	mid = (a + b) / 2
	err = (b - a) / 2

	fmid = func(mid)
	fa = func(a)
	fb = func(b)

	while abs(fmid) > tol

		# println("a = $a     mid = $mid     b = $b")
		# println("fa = $fa     fmid = $fmid     fb = $fb")

		if fa * fmid < 0
			b = mid
			fb = func(b)

		elseif fmid * fb < 0
			a = mid
			fa = func(a)

		elseif fmid == 0.0
			err = 0
			return (mid, err)

		end

		mid = (a + b) / 2
		err = (b - a) / 2
		fmid = func(mid)

		if mid == a || mid == b
			return (mid, err)
		end

	end

	return (mid, err)
end

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
export simp
# Simpson method integration
function simp()

end


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
export countin
# Return dict with number of repetitionsof elements in the array
function countin(arr, acc=16)
	dct = Dict{Float64, Int64}()
	for i = 1:length(arr)
		key = round(arr[i], acc);
		if haskey(dct, key)
			dct[key] += 1;
		else 
			dct[key] = 1;
		end
		# Workaround on -0.0 and +0.0
		if haskey(dct, -0.0)
			if haskey(dct, 0.0)
				dct[0.0] += dct[-0.0]
			else 
				dct[0.0] = dct[-0.0]
			end
			delete!(dct, -0.0)
		end
	end
	return dct;	
end

export countin2
function countin2(arr, acc=16)
	dct = countin(arr, acc)
	x = sort( [k for k in keys(dct)] )
	y = Array(Int64, length(dct))
	for i = 1:length(x)
		y[i] = dct[x[i]];
	end
	return x, y;
end

end
