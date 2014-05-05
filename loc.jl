#!/usr/bin/env julia

module loc


using Distributions



type Site 
	coord 		::Array{Float64,1}
	energy		::Float64
	checkouts 	::Int32
end



function hypot3(x, y, z)
	return hypot(hypot(x, y), z)
end


# Distance between sites as dots for three coordinates
function distanceXYZ(s1::Site, s2::Site)
	return hypot3(	s1.x - s2.x,
					s1.y - s2.y,
					s1.z - s2.z  )
end




# Recursive assignments of energy and calculation of ...
function walk_and_assign(domain, i, j, k, step, 
						checkouts )

end




# Gather the statistics on recombination of excitons 
# in normally distributed potential
# gather(	iter_num, pot_form_num, domain_size, grid_step, distr )
function gather(	iter_num	::Int32, 
					pot_form_num::Int32,
					domain_size	::Array{Int32,1},
					grid_step	::Float64
					distr		::Distribution )
	
	# domain = Array(Site, domain_size, domain_size, domain_size)
	domain = ...
	
	# Iterate over sites
	# Generate all sites in the field or do it on-the-fly?
	for domain_index = 1:pot_form_num

		# Iterate over the given/generating sites array
		for iteration_index = 1:iter_num

			# Initial indexes
			(i, j, k) = int(domain_size * rand(3)) 
			

			# Work on recursion on sites

		end
	end
end






end