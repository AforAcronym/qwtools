#!/usr/bin/env julia

module loc

export gather
export ev_div_kt
export erg_div_kt

using utils
using constants
using Distributions



type Site 
    pos    ::Array{Float64,1}
    energy ::Float64
end



# # function hypot3(x, y, z)
# #     return hypot(hypot(x, y), z)
# # end


# # # Distance between sites as dots for three coordinates
# # function distanceXYZ(s1::Site, s2::Site)
# #     return hypot3(  s1.x - s2.x,
# #                     s1.y - s2.y,
# #                     s1.z - s2.z  )
# # end





# # function populate_perimeter(arr, distr::Distribution)
# #     (rows, cols, depz) = size(arr)
# #     arr[ 2:end,   1,       :    ] = rand( d, length( arr[2:end,   1,       :    ]) )
# #     arr[   1,     :,       :    ] = rand( d, length( arr[  1,     :,       :    ]) )
# #     arr[ 2:end, 2:end,     1    ] = rand( d, length( arr[2:end, 2:end,     1    ]) )
# #     arr[ 2:end, 2:end,    end   ] = rand( d, length( arr[2:end, 2:end,    end   ]) )
# #     arr[ 2:end,  end,    2:end-1] = rand( d, length( arr[2:end,  end,    2:end-1]) )
# #     arr[  end,  2:end-1, 2:end-1] = rand( d, length( arr[ end,  2:end-1, 2:end-1]) )
# # end


# # # Requires function of 1 argument
# # function populate_perimeter(arr, f::Function)
# #     arr[ 2:end,   1,       :    ] = [f(x) for x in arr[2:end,   1,       :    ]]
# #     arr[   1,     :,       :    ] = [f(x) for x in arr[  1,     :,       :    ]]
# #     arr[ 2:end, 2:end,     1    ] = [f(x) for x in arr[2:end, 2:end,     1    ]]
# #     arr[ 2:end, 2:end,    end   ] = [f(x) for x in arr[2:end, 2:end,    end   ]]
# #     arr[ 2:end,  end,    2:end-1] = [f(x) for x in arr[2:end,  end,    2:end-1]]
# #     arr[  end,  2:end-1, 2:end-1] = [f(x) for x in arr[ end,  2:end-1, 2:end-1]]
# # end


# # function populate_perimeter(arr::arr{Float64,3}, 
# #                             distr::Distribution, 
# #                             exc_lifetime::Float64)
# #     (rows, cols, depz) = size(arr)
# #     arr[ 2:end,   1,       :    ] = rand( d, length( arr[2:end,   1,       :    ]) )
# #     arr[   1,     :,       :    ] = rand( d, length( arr[  1,     :,       :    ]) )
# #     arr[ 2:end, 2:end,     1    ] = rand( d, length( arr[2:end, 2:end,     1    ]) )
# #     arr[ 2:end, 2:end,    end   ] = rand( d, length( arr[2:end, 2:end,    end   ]) )
# #     arr[ 2:end,  end,    2:end-1] = rand( d, length( arr[2:end,  end,    2:end-1]) )
# #     arr[  end,  2:end-1, 2:end-1] = rand( d, length( arr[ end,  2:end-1, 2:end-1]) )
# # end








# Position shift needed to reduce indexes
function ipos_shift( x, xlim )
    x < 1    && return  xlim + ipos_shift(x + xlim, xlim)
    x > xlim && return -xlim + ipos_shift(x - xlim, xlim)
    return 0
end






# Relative decay
# domain to calculate
# ipos — index position of the current site in the domain
# lim  — maximum number of sites layers around to search  within
function rel_decay( domain::Array{Site, 3}, ipos::Array{Int64,1}, lim::Int64 )
        
        step = domain[1,1,1].pos[1] - domain[2,1,1].pos[1]
        
        (rows, cols, depz) = size(domain)

        # Dimensions length for the desired cube
        dim_len = 2 * lim + 1 
        rel_decays = Array(Float64, dim_len, dim_len, dim_len)
        
        shift = ipos - [1,1,1]   # Indexes shift from the origin

        for d = ipos[3]-lim:ipos[3]+lim, c = ipos[2]-lim:ipos[2]+lim, r = ipos[1]-lim:ipos[1]+lim

            # println("ipos = ", ipos)
            # println("shift = ", shift)
            # println("lim = ", lim)
            # println("[r, c, d] = ", [r, c, d])

            # Shifts for subindexes
            rshft = ipos_shift(r, rows)
            cshft = ipos_shift(c, cols)
            dshft = ipos_shift(d, depz)

            # println("[r+rshft, c+cshft, d+dshft] = ", [r+rshft, c+cshft, d+dshft])

            dest_site = domain[r + rshft, c + cshft, d + dshft]
            orig_site = domain[ipos...] # Check

            # distance_ratio = vecnorm(dest_site.pos - [rshft, cshft, dshft].*step - orig_site.pos)
            distance_ratio = vecnorm([r, c, d].*step - orig_site.pos)
            # println("dest: ", dest_site.pos)
            # println("shift: ", [rshft, cshft, dshft].*step)
            # println("orig: ", orig_site.pos)
            # println("distance: ", distance_ratio)

            energy_ratio = 0 
            if dest_site.energy > orig_site.energy
                energy_ratio = (dest_site.energy - orig_site.energy)
            end
            # println("[r,c,d] - shift .+ lim = ", [r,c,d] - shift .+ lim)
            
            rel_decays[ ([r,c,d] - shift .+ lim)... ] = exp( -distance_ratio - energy_ratio ) # Check
            
            
            if isnan(exp( -distance_ratio - energy_ratio ))
                println("-distance_ratio = ", -distance_ratio)
                println("-energy_ratio = ", -energy_ratio)
                println("-distance_ratio - energy_ratio  = ", -distance_ratio - energy_ratio )
                println("ipos = ", ipos)
                println("shift = ", shift)
                println("lim = ", lim)
                println("[r, c, d] = ", [r, c, d])
                println("[r+rshft, c+cshft, d+dshft] = ", [r+rshft, c+cshft, d+dshft])
                println("dest: ", dest_site.pos)
                println("shift: ", [rshft, cshft, dshft].*step)
                println("orig: ", orig_site.pos)
                println("distance: ", distance_ratio)
                println("[r,c,d] - shift .+ lim = ", [r,c,d] - shift .+ lim)
            end

        end

        # rel_decays[lim, lim, lim] = 0
        return rel_decays


end




# 

# total_time, Site = jump(domain, [i, j, k], exc_lifetime, temperature, total_time, lim)
function jump(  domain      ::Array{Site, 3}, 
                ipos        ::Array{Int64,1}, 
                exc_lifetime::Float64, 
                total_time  ::Float64,
                lim         ::Int64  )

    # Taking step into account, estimate number of envelope perimeters
    # sufficient for the calculation. No need to go through sites which stand
    # far enough to be neglected. lim is now used for the purpose

    # This decays are relative because they are not multiplied with
    # exc_lifetime since it is required only for the hop_time
    rel_decays = rel_decay(domain, ipos, lim)
    println("Sum of rel_decays: ", sum(rel_decays) )
    rel_decay_total = sum(rel_decays) + 1.0

    hop_time = - exc_lifetime / log(rand()) / rel_decay_total   # WTF? Baranovskii PysRevB 58, 19 (1998)
    # hop_time = - exc_lifetime * log(rand()) / rel_decay_total # WTF? Schoenherr ChemPhys 52, 287 (1980)

    # Hop
    if hop_time < exc_lifetime 
        
        # probs = rel_decays ./ rel_decay_total               # Probabilities of hops, sum(probs) < 1
        probs = rel_decays ./ sum(rel_decays)                 # Probabilities of hops — sum = 1
        println("probs sum: ", sum(probs))
        probs_distribution = Categorical(vec(probs))          # Probabilities distribution
        
        next_site_of_choise = rand( probs_distribution )      # Choose a site to jump to (index)
        next_ipos = ind2sub(size(probs), next_site_of_choise) # Obtain (i, j, k) from the index
        # println("next_ipos: ", next_site_of_choise, " ", next_ipos)

        # Current site is chosen to be the next XXX check
        # Recombination
        if ipos == next_ipos
            # println("wow, staying here") 
            return total_time + exc_lifetime, domain[ipos...]
        end

        return jump(domain, [next_ipos...], exc_lifetime, total_time + hop_time, lim)
    end

    # Recombination
    return total_time + exc_lifetime, domain[ipos...]
end






# Calculate E/kT, E in electron-volts, T in kelvins
function ev_div_kt(energy, temp)
    return energy / convert_erg2ev(CONST_BOLTZMANN) / temp
end




# Calculate E/kT, E in ergs, T in kelvins
function erg_div_kt(energy, temp)
    return energy / CONST_BOLTZMANN / temp
end






# Gather the statistics on recombination of excitons in arbitrary distributed
# potential. Since real physical units don't play role here, whereas only 
# their ratios do, follow this:
#   - grid_step should be lattice parameter normed by localization length 
#     (lat_par/loc_len)
#   - pot_distr should be normed by thermal energy (E/kT)
function gather( iter_num     ::Int64, 
                 pot_form_num ::Int64,
                 domain_size  ::Array{Int64,1},
                 grid_step    ::Float64,
                 pot_distr    ::Distribution,
                 exc_lifetime ::Float64,
                 lim          ::Int64     )
    
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
        domain[r, c, d] = Site([r, c, d] .* grid_step, 0)
    end

    # Iterate over the given/generating sites array
    for index_dom = 1:pot_form_num       
        
        # Populate energy
        for i = 1:length(domain)
            domain[i].energy = rand(pot_distr)
        end
        
        # Place excitons for hopping
        for index_exc = 1:iter_num         

            # Random initial indexes, avoid zeros
            i, j, k = rand(1:rows), rand(1:cols), rand(1:depz)

            # Start the hopping recursion over the domain
            decay_time::Float64, last_site::Site = jump(domain, [i,j,k], exc_lifetime, 0.0, lim)

            decay_times[index_exc, index_dom] = decay_time
            last_sites[ index_exc, index_dom] = last_site
        end
    end

    return decay_times, last_sites
end






end