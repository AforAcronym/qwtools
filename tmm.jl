#!/usr/bin/env julia


module tmm

export Region
export Heterostructure
export test
export eigenenergy


using constants
using utils
# using Winston



const CALC_TOLERANCE    = eps(Float64)
const CALC_STEP         = convert_ev2erg(0.001)



type Region
    potential::Float64
    width::Float64
    center::Float64
    mass::Float64

    Region(potential, width, center, mass) = 
        width < 0 ? error("Region cannot have negative width") : new(potential, width, center, mass)
    Region(potential, width, center, mass) =  
        mass < 0 ? error("Region cannot have negative mass")  : new(potential, width, center, mass)
end



immutable Heterostructure
    layers::Array{Region}
end








# Index stands here for the left region and matrix is calculated for the propagation from 
# this region to the next one.
function propagationMatrix(hs::Heterostructure, index::Int, energy::Float64)

    if index >= length(hs.layers)
        println("index = $index, heterostructure length = $length(hs.layers)")
        println("Energy = $energy")
        error("Index requested for propagation matrix is out of bounds")
    end

    # Propagation from index (1) to index+1 (2)
    m1sqrt   = sqrt( hs.layers[ index ].mass )
    m2sqrt   = sqrt( hs.layers[index+1].mass )
    EmV1sqrt = sqrt( complex(energy - hs.layers[ index ].potential) )
    EmV2sqrt = sqrt( complex(energy - hs.layers[index+1].potential) )
    w1       = hs.layers[ index ].width
    w2       = hs.layers[index+1].width

    ampl     = (m2sqrt / m1sqrt) * (EmV1sqrt / EmV2sqrt)
    ampl_pls = 0.5 * (1 + ampl)
    ampl_mns = 0.5 * (1 - ampl)

    phase_coef = im * sqrt(ELECTRON_MASS / 2) / CONST_PLANCK_REDUCED
    phase1     = w1 * m1sqrt * EmV1sqrt
    phase2     = w2 * m2sqrt * EmV2sqrt
    phase_pls  = phase_coef * (phase2 + phase1)
    phase_mns  = phase_coef * (phase2 - phase1)

    return [ ampl_pls * exp( phase_pls)     ampl_mns * exp( phase_mns)
             ampl_mns * exp(-phase_mns)     ampl_pls * exp(-phase_pls) ]

end






function propagationAmpl(hs::Heterostructure, energy::Float64)
    m = propagationMatrix(hs, 1, energy)
    for i = 2 : length(hs.layers) - 1
        m = propagationMatrix(hs, i, energy) * m
    end
    # TODO calculate and add A and B amplitudes to hs for wavefunctions building?
    # Save an array of amplitudes A and B. If the calculated energy is eigen, output it?
    # ...or calculate the amplitudes separately?
    return abs2(m[2,2]) # abs2() is faster than abs()
end






function propagAmplDeriv(hs::Heterostructure, energy::Float64)
    step = eps(Float64)
    return (log1p(propagationAmpl(hs, energy + step)) - 
            log1p(propagationAmpl(hs, energy - step))) / step / 2
end




function energyRangeInfo(hs::Heterostructure)
    len = size(hs.layers, 1)
    potentials = Array(Float64, len)
    for i = 1:len
        potentials[i] = hs.layers[i].potential 
    end
    offset = minimum(potentials)
    top = sort(potentials)[end-1]
    return (top, offset)
end




function eigenenergy(hs::Heterostructure, energy_step::Float64, number_of_states::Int=1000)

    if number_of_states <= 0
        error("There is no point in calculation for zero or negative number of states")
    end 
    
    energy_eigen = Float64[] #Array(Float64, 1)
    energy_error = Float64[] #Array(Float64, 1)
    energy_max, energy = energyRangeInfo(hs)
    energy = energy + energy_step
    energy_next = energy + energy_step

    zerofunc(x) = propagAmplDeriv(hs, x)

    while energy < energy_max
        if energy_next > energy_max; energy_next = energy_max; end
        
        if ((propagAmplDeriv(hs, energy) * propagAmplDeriv(hs, energy + energy_step) < 0) 
            && propagAmplDeriv(hs, energy) < 0 )
        
            (energy_extremum, err) = bisection(zerofunc, energy, energy + energy_step, CALC_TOLERANCE)
            
            push!(energy_eigen, energy_extremum) 
            push!(energy_error, err)
        end

        if length(energy_eigen) == number_of_states; break; end

        energy = energy_next
        energy_next = energy_next + energy_step

    end

    return energy_eigen, energy_error
end


# 
function eigenenergy(hs::Heterostructure)
    eigenenergy(hs::Heterostructure, CALC_STEP) 
end



function eigenenergy(hs::Heterostructure, step::Float64)
    eigenenergy(hs, step, 1000)
end



function eigenenergy(hs::Heterostructure, number_of_states::Int=1000)
    eigenenergy(hs, CALC_STEP, number_of_states)
end




function test(roundto=3)
    p = convert_ev2erg(0.4)     # potential, eV -> erg
    # ww = 3e-7                   # quantum well width, 3 nm
    ww = 3e-6                   # quantum well width, 3 nm
    wb = 2e-5                   # barrier width, 20 nm
    m = 0.067                   # particle mass

    hs = Heterostructure([  Region(p, wb, -(ww+wb)/2, m) 
                            Region(0, ww,   0,        m)
                            Region(p, wb,  (ww+wb)/2, m) ])
    tic()
    energies, errors = eigenenergy(hs, convert_ev2erg(0.001), 9) # Resolution 1 meV
    toc()

    energies_ev = convert_erg2ev(energies)
    # errors_ev   = convert_erg2ev(errors)

    # println( join( ["$i: $(energies_ev[i]) ± $(errors_ev[i]) eV" for i = 1:length(energies_ev)] ,"\n") )
    println( join( ["$i: $(round(energies_ev[i],roundto)) eV" for i = 1:length(energies_ev)] ,"\n") )

end


end
