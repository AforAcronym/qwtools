#!/usr/bin/env julia

module loc

export gather
export ev_div_kt
export erg_div_kt
export test_stats

using utils
using constants
using Distributions
using PyPlot


type Site 
    pos    ::Array{Float64,1}
    energy ::Float64
end



# Position shift needed to reduce indexes
function ipos_shift( x, xlim )
    x < 1    && return  xlim + ipos_shift(x + xlim, xlim)
    x > xlim && return -xlim + ipos_shift(x - xlim, xlim)
    return 0
end





# Dimensions length for the desired decays cube
function dim_length(domsize::Array{Int64,1}, lim::Int64)
    if length(domsize) != 3
        error("Domain must be three-dimensional.")
    end
    dim_len = [lim,lim,lim] .* 2 .+ 1
    for i = 1:3
        # Choose minimum for dimension length
        dim_len[i] = domsize[i] < dim_len[i] ? domsize[i]      : dim_len[i]
        # dim_len must be odd (even-1) in order to be symmetric
        dim_len[i] = dim_len[i] % 2 == 0     ? dim_len[i] - 1  : dim_len[i]
    end
    lims = int((dim_len .- 1) ./ 2)
    return dim_len, lims
end





# Relative decay
# domain to calculate
# ipos — index position of the current site in the domain
# lim  — maximum number of sites layers around to search  within
function rel_decay( domain::Array{Site, 3}, ipos::Array{Int64,1}, lim::Int64 )
        # println("-----------------------")
        # println("-----------------------")
        # println("-----------------------")
        
        step = domain[1].pos[1] - domain[2].pos[1]
        (rows, cols, depz) = size(domain)

        # Dimensions length for the desired cube
        dim_len, lims = dim_length([rows, cols, depz], lim)

        rel_decays = Array(Float64, dim_len...)
        
        shift = ipos - [1,1,1]   # Indexes shift from the origin
        

        for d = ipos[3] - lims[3] : ipos[3] + lims[3], 
            c = ipos[2] - lims[2] : ipos[2] + lims[2], 
            r = ipos[1] - lims[1] : ipos[1] + lims[1]

            # println("- - - - - - - - ")
            # println("ipos    = ", ipos)
            # println("shift   = ", shift)
            # println("lim     = ", lim)
            # println("[r,c,d] = ", [r, c, d])

            # Shifts for subindexes
            rshft = ipos_shift(r, rows)
            cshft = ipos_shift(c, cols)
            dshft = ipos_shift(d, depz)

            # println("[r+rshft, c+cshft, d+dshft] = ", [r+rshft, c+cshft, d+dshft])

            dest_site = domain[r + rshft, c + cshft, d + dshft]
            orig_site = domain[ipos...]

            # distance_ratio = vecnorm(dest_site.pos - [rshft, cshft, dshft].*step - orig_site.pos)
            distance_ratio = vecnorm([r, c, d].*step - orig_site.pos)
            # println("dest: ", dest_site.pos)
            # println("ishift: ", [rshft, cshft, dshft].*step)
            # println("orig: ", orig_site.pos)
            # println("distance: ", distance_ratio)

            energy_ratio = 0 
            if dest_site.energy > orig_site.energy
                energy_ratio = (dest_site.energy - orig_site.energy)
            end
            # println("dim_len: ", dim_len)
            # println("lims: ", lims)
            # println("[r,c,d] - shift .+ lim = ", [r,c,d] - shift .+ lim)
            # println("[r + lims[1], c + lims[2], d + lims[3]] - shift = ", [r + lims[1], c + lims[2], d + lims[3]] - shift)
            
            rel_decays[ ([r + lims[1], c + lims[2], d + lims[3]] - shift)... ] = exp( -distance_ratio - energy_ratio )
            # rel_decays[ ([r, c, d] - shift .+ lim)... ] = exp( -distance_ratio - energy_ratio )

        end

        # rel_decays[lim, lim, lim] = 0
        return rel_decays


end




# 

# total_time, Site = jump(domain, [i, j, k], exc_lifetime, temperature, total_time, lim)
function jump(  domain      ::Array{Site, 3}, 
                ipos        ::Array{Int64,1}, 
                exc_lifetime::Float64, 
                esc_rate    ::Float64, 
                total_time  ::Float64,
                lim         ::Int64  )

    # This decays are relative because they are not multiplied by exc_lifetime
    decays = rel_decay(domain, ipos, lim) .* esc_rate # v(i,j)
    
    # println("Sum of rel_decays: ", sum(rel_decays) )
    # println("rel_decays[lim,lim,lim]: ", rel_decays[lim,lim,lim] )
    
    _, lims = dim_length([size(domain)...], lim)
    center = lims .+ 1
    # println("lims: ", lims)
    decay_current = sum(decays) - decays[center...] + (exc_lifetime)^-1 # v(i)
    
    # # # # # # # # # # !!!!!!!!!!!!!!!!!!!!
    # decays[center...] = exc_lifetime^-1 
    decays[center...] = 0 

    hop_time = -1 / log(rand()) / decay_current   # t(i), Baranovskii PysRevB 58, 19 (1998)
    # XXX NOTE WTF? Schoenherr ChemPhys 52, 287 (1980): -1 * log(rand()) / decay_current

    # Hop
    if hop_time < exc_lifetime 
        
        probs = decays ./ decay_current               # Probabilities of hops, P(i,j), sum(probs) < 1
        probs ./= sum(probs)                          # Technically correct normalization
        # println("decays sum: ", sum(decays))
        # println("decays*esc_rate sum: ", sum(decays) * esc_rate)
        # println("probs sum: ", sum(probs))
        # probs[center...] = 1.0 - sum(probs)         # FIXME Why?
        
        probs_distribution = Categorical(vec(probs) ./ sum(probs))  # Probabilities distribution
            
        next_site =  rand( probs_distribution )       # Choose a site to jump to (index)
        next_ipos = ind2sub(size(probs), next_site)   # Obtain (i, j, k) from the index
        # println("next_ipos: ", next_site, " ", next_ipos)

        # Current site is chosen to be the next XXX check
        # Recombination
        if ipos == next_ipos
            println("wow, staying here") 
            return total_time + exc_lifetime, domain[ipos...]
        end

        return jump(domain, [next_ipos...], exc_lifetime, esc_rate, total_time + hop_time, lim)
    end

    # Recombination
    return total_time + exc_lifetime, domain[ipos...]
end






# Calculate E/kT, E in electron-volts, T in kelvins
function ev_div_kt(energy_ev, temperature)
    return energy_ev / convert_erg2ev(CONST_BOLTZMANN) / temperature
end




# Calculate E/kT, E in ergs, T in kelvins
function erg_div_kt(energy_erg, temperature)
    return energy_erg / CONST_BOLTZMANN / temperature
end






# Gather the statistics on recombination of excitons in arbitrary distributed
# potential. Since real physical units don't play role here, whereas only 
# their ratios do, follow this:
#   - grid_step should be lattice parameter normed by localization length 
#     (lat_par/loc_len)
#   - pot_distr should be normed by thermal energy (E/kT)
# gather( 10, 10, [70,70,70], 1.0, Normal(0, ev_div_kt(0.05, 5)), 1e-13, 1e16, 6)
function gather( iter_num     ::Int64, 
                 pot_form_num ::Int64,
                 domain_size  ::Array{Int64,1},
                 grid_step    ::Float64,
                 pot_distr    ::Distribution,
                 exc_lifetime ::Float64,
                 esc_rate     ::Float64,
                 lim          ::Int64     )
    
    if length(domain_size) != 3
        error("Domain must be 3-dimensional, e.g. [20 40 50] or [100 200 1] in case of 2D")
    end
    
    (rows, cols, depz) = domain_size

    # Allocation
    domain      = Array(Site,    rows, cols, depz)
    # decay_times = Array(Float64, iter_num, pot_form_num)
    # exit_sites  = Array(Site,    iter_num, pot_form_num)
    decay_times = Dict{Float64, Int64}()
    exit_energy  = Dict{Float64, Int64}()

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
            decay_time::Float64, exit_site::Site = jump(domain, [i,j,k], exc_lifetime, esc_rate, 0.0, lim)

            # decay_times[index_exc, index_dom] = decay_time
            # exit_sites[ index_exc, index_dom] = exit_site
            add2dict!(decay_times, decay_time)
            add2dict!(exit_energy, exit_site.energy)
        end
    end

    return decay_times, exit_energy
end





# include("loc.jl"); @time stats = loc.test_stats(1.0, 300.0, 74.0, 
#                    40, 30, [200,200,200], 
#                    10.0, 0.05, 1e-13, 1e16, 6); loc.stats_plot(stats)
# 
# include("loc.jl"); @time stats = loc.collect_stats(10, 10, [70,70,70], 0.5, 
# loc.ev_div_kt(0.05, 5), 1e-13, 1e16, 6, float64([5,20,50,100,200,300])); loc.plot_stats(stats, 3)
function collect_stats(iter_num::Int64             = 10,
                       pot_form_num::Int64         = 5,
                       domain_size::Array{Int64,1} = [70,70,70],
                       grid_step::Float64          = 5.0,
                       sigma::Float64              = 0.05,
                       exc_lifetime::Float64       = 1e-13,
                       esc_rate::Float64           = 1e16,
                       lim::Int64                  = 6,
                       temps::Array{Float64,1}     = float64([  
                                                    1,   5,  10,  15,  20,  25,  30,  40,  
                                                   50,  60,  70,  80,  90, 100, 110, 125, 
                                                  140, 155, 170, 200, 230, 260, 290 ])      
                       )
    
    stats = Dict{ Float64, Dict{String, Dict{Float64, Int64}} }()
    counter = 0
    countmax = length(temps)
    for t in temps
        
        counter += 1
        println(counter, "/", countmax, ": Calculation for T = ", round(t,1))

        @time decay_times, exit_energy = gather(iter_num,    pot_form_num,
                                               domain_size,  grid_step,
                                               Normal(0, ev_div_kt(sigma, t)), 
                                               exc_lifetime, esc_rate, lim)
        stats[t] = [ "decay_times" => decay_times, 
                     "exit_energy" => exit_energy  ]
    end
    
    return stats
        
end





# Electron-volts!
# FIXME Make peaks instead of averaging lines?
function plot_stats( stats, roundto::Int64=3 )
    fig = figure("Simulated PL")
    for t in sort([k for k in keys(stats)])
        if t == 0; continue; end # bad
        # decay_times = stats[t]["decay_time"] # TODO
        # exit_energy  = stats[t]["exit_energy"]
        energy = [ nrg*convert_erg2ev(CONST_BOLTZMANN)*t 
                    for nrg in keys(stats[t]["exit_energy"])]
        # counts_dict = {nrg => 0 for nrg in unique(energy)}
        # for nrg in energy
        #     counts_dict[nrg] += 1
        # end
        # energy_sorted = sort([k for k in keys(counts_dict)])
        # counts = [counts_dict[nrg] for nrg in energy_sorted]
        # len = length(counts)
        # plot3D(ones(len) .* t, energy_sorted, counts)
        x, y = countin2(energy, roundto)
        plot3D(ones(length(x)) .* t, x, y)
    end
end



end