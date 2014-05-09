#!/usr/bin/env julia

module loc

using utils
using Distributions



type Site 
    # coord    ::Vector{Float64}
    coord    ::Array{Float64,1}
    energy   ::Float64
end



# function hypot3(x, y, z)
#     return hypot(hypot(x, y), z)
# end


# # Distance between sites as dots for three coordinates
# function distanceXYZ(s1::Site, s2::Site)
#     return hypot3(  s1.x - s2.x,
#                     s1.y - s2.y,
#                     s1.z - s2.z  )
# end





function populate_perimeter(arr, distr::Distribution)
    (rows, cols, depz) = size(arr)
    arr[ 2:end,   1,       :    ] = rand( d, length( arr[2:end,   1,       :    ]) )
    arr[   1,     :,       :    ] = rand( d, length( arr[  1,     :,       :    ]) )
    arr[ 2:end, 2:end,     1    ] = rand( d, length( arr[2:end, 2:end,     1    ]) )
    arr[ 2:end, 2:end,    end   ] = rand( d, length( arr[2:end, 2:end,    end   ]) )
    arr[ 2:end,  end,    2:end-1] = rand( d, length( arr[2:end,  end,    2:end-1]) )
    arr[  end,  2:end-1, 2:end-1] = rand( d, length( arr[ end,  2:end-1, 2:end-1]) )
end


# Requires function of 1 argument
function populate_perimeter(arr, f::Function)
    arr[ 2:end,   1,       :    ] = [f(x) for x in arr[2:end,   1,       :    ]]
    arr[   1,     :,       :    ] = [f(x) for x in arr[  1,     :,       :    ]]
    arr[ 2:end, 2:end,     1    ] = [f(x) for x in arr[2:end, 2:end,     1    ]]
    arr[ 2:end, 2:end,    end   ] = [f(x) for x in arr[2:end, 2:end,    end   ]]
    arr[ 2:end,  end,    2:end-1] = [f(x) for x in arr[2:end,  end,    2:end-1]]
    arr[  end,  2:end-1, 2:end-1] = [f(x) for x in arr[ end,  2:end-1, 2:end-1]]
end


function populate_perimeter(arr::arr{Float64,3}, 
                            distr::Distribution, 
                            exc_lifetime::Float64)
    (rows, cols, depz) = size(arr)
    arr[ 2:end,   1,       :    ] = rand( d, length( arr[2:end,   1,       :    ]) )
    arr[   1,     :,       :    ] = rand( d, length( arr[  1,     :,       :    ]) )
    arr[ 2:end, 2:end,     1    ] = rand( d, length( arr[2:end, 2:end,     1    ]) )
    arr[ 2:end, 2:end,    end   ] = rand( d, length( arr[2:end, 2:end,    end   ]) )
    arr[ 2:end,  end,    2:end-1] = rand( d, length( arr[2:end,  end,    2:end-1]) )
    arr[  end,  2:end-1, 2:end-1] = rand( d, length( arr[ end,  2:end-1, 2:end-1]) )
end







# Transition rate between sites
# Energy is in eV
function rate(s1::Site, s2::Site, exc_lifetime, temperature)
    energy_ratio = 0
    if s2.energy > s1.energy  
        energy_ratio = (s2.energy - s1.energy) / temperature /  convert_erg2ev(CONST_BOLTZMANN)
    end
    # Baranovskii
    return exp( -vecnorm(s1.coord - s2.coord) - energy_ratio) / exc_lifetime
end


function pos_shift( x, xlim )
    if x < 0
        return xlim
    elseif x > xlim
        return -xlim
    end
    return 0
end







function rel_decay( domain     ::Array{Site, 3},
                    ipos       ::Array{Int64,1},
                    temperature::Float64  )
        
        step = domain[1,1,1].coord[1] - domain[2,1,1].coord[1]
        
        (rows, cols, depz) = size(domain)

        # Limit for number of lattice layers to research
        lim = 4                 # FIXME magic number
        
        dim_len = 2 * lim + 1
        rel_decays = Array(Float64, dim_len, dim_len, dim_len)
        
        shift = ipos - [1,1,1]   # Shift from the origin

        for d = ipos[1]-lim:ipos[1]+lim, c = ipos[2]-lim:ipos[2]+lim, r = ipos[3]-lim:ipos[3]+lim

            # Shifts for subindexes
            rshft = ipos_shift(r, rows)
            cshft = ipos_shift(c, cols)
            dshft = ipos_shift(d, depz)

            dest_site = domain[r + rshft, c + cshft, d + dshft]
            orig_site = domain[ipos...]

            distance_ratio = vecnorm(dest_site.coord - [rshft, cshft, dshft].*step - orig_site.coord)  
            
            energy_ratio = 0 
            if dest_site.energy > orig_site.energy  
                energy_ratio = (dest_site.energy - orig_site.energy) / 
                                                (temperature * convert_erg2ev(CONST_BOLTZMANN))
            end

            # FIXME Incosistent approach to distance (no unit) and energy (has unit)
            rel_decays[ ([r,c,d] - shift .+ lim)... ] = exp( -distance_ratio - energy_ratio ) # Check
        end

        rel_decays[lim, lim, lim] = 0
        return rel_decays


end






# total_time, Site = jump(domain, [i, j, k], exc_lifetime, temperature, total_time)
function jump(  domain      ::Array{Site, 3}, 
                pos         ::Array{Int64,1}, 
                exc_lifetime::Float64, 
                temperature ::Float64,
                total_time  ::Float64  )

    # Taking step into account, estimate number of envelope perimeters
    # sufficient for the calculation. No need to go through sites
    # which stand far enough to be neglected

    # This decays are relative because they are not multiplied with
    # exc_lifetime since it is required only for the hop_time
    rel_decays = rel_decay(domain, pos, temperature) # not implemented
    rel_decay_total = sum(rel_decays) + 1.0

    hop_time = -log(rand()) / rel_decay_total * exc_lifetime

    # Hop
    if hop_time < exc_lifetime 
        probs = rel_decays ./ rel_decay_total                # Probabilities of hops
        probs_distribution = Categorial(vec(probs))          # Probabilities distribution
        next_site_of_choise = rand( prob_distribution )      # Choose a site to jump to (index)
        next_pos = ind2sub(size(probs), next_site_of_choise) # Obtain [i, j, k] from the index
        return jump(domain, next_pos, exc_lifetime, temperature, total_time + hop_time)
    end

    # Recombination
    return total_time + exc_lifetime, domain[pos...]
end






# Gather the statistics on recombination of excitons in arbitrary distributed
# potential
function gather( iter_num     ::Int32, 
                 pot_form_num ::Int32,
                 domain_size  ::Array{Int32,1},
                 grid_step    ::Float64
                 pot_distr    ::Distribution,
                 exc_lifetime ::Float64,
                 temperature  ::Float64,
                 loc_length   ::Float64   )
    
    if length(domain_size) != 3
        error("Domain must be 3-dimensional, e.g. [20 40 50] or [20 40 1] in case of 2D")
    end
    
    (rows, cols, depz) = domain_size

    # Allocation
    domain      = Array(Site,    rows, cols, depz)
    decay_times = Array(Float64, iter_num, pot_form_num)
    last_sites  = Array(Site,    iter_num, pot_form_num)

    # Populate coordinates
    for d = 1:depz, c = 1:cols, r = 1:rows
        # Take scaling factor into account
        domain[r, c, d] = Site([r c d] .* (grid_step / loc_length), 0)
    end

    # Iterate over the given/generating sites array
    for index_dom = 1:pot_form_num       
        
        # Populate energy
        for i = 1:length(domain)
            domain[i].energy = rand(pot_distr)
        end
        
        # Place excitons for hopping
        for index_exc = 1:iter_num         

            # Random initial indexes
            (i, j, k) = int(domain_size * rand(size(domain_size))) 

            # Start the hopping recursion over the domain
            decay_time::Float64, last_site::Site = jump(domain, [i,j,k], exc_lifetime, temperature, 0)

            # Take scaling factor into account again
            last_site.coord .*= loc_length

            decay_times[index_exc, index_dom] = decay_time
            last_sites[ index_exc, index_dom] = last_site
        end
    end

    return decay_times, last_sites
end






end