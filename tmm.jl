#!/usr/bin/env julia


module tmm

export Region
export Heterostructure
export test
export eigenenergy
export waveFunction
export drawwf

using constants
using utils
using Roots

using PyPlot
using PyCall
@pyimport numpy



const ERG_CALC_STEP         = convert_ev2erg(0.001)
const CALC_TOLERANCE        = eps(Float64)
# const CALC_TOLERANCE  = ERG_CALC_STEP / 100

# FIXME fet rid of it, see drawwf() function below
const WAVE_FUNCTION_DIVIDER = 2e16 



type Region
    potential::Float64
    width    ::Float64
    center   ::Float64
    mass     ::Float64

    Region(potential, width, center, mass) = 
        width < 0 ? 
            error("Region cannot have negative width") : 
            new(potential, width, center, mass)

    Region(potential, width, center, mass) =  
        mass < 0 ? 
            error("Region cannot have negative mass") : 
            new(potential, width, center, mass)
end



immutable Heterostructure
    layers::Array{Region}
    energy_levels::Array{Float64}
end





function waveNumber(layer::Region, energy::Float64)
    return sqrt(complex(energy - layer.potential) * layer.mass * ELECTRON_MASS * 2) / 
                                                            CONST_PLANCK_REDUCED
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





# Set of all the transfer matrices, can be unused
function propagationMatrixSet(hs::Heterostructure, energy::Float64)
    return [propagationMatrix(hs, i, energy) for i in length(hs.layers)-1:-1:1]
end





function propagationAmpl(hs::Heterostructure, energy::Float64)
    m = propagationMatrix(hs, 1, energy)
    for i = 2 : length(hs.layers) - 1
        m = propagationMatrix(hs, i, energy) * m
    end
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




function eigenenergy(   hs              ::Heterostructure, 
                        energy_step     ::Float64, 
                        number_of_states::Int=1000  )

    if number_of_states <= 0
        error("There is no point in calculation for zero or negative number of states")
    end 

    energy_max, energy = energyRangeInfo(hs)
    energy = energy + energy_step
    energy_next = energy + energy_step

    zerofunc(x) = propagAmplDeriv(hs, x)

    while energy < energy_max
        if energy_next > energy_max; energy_next = energy_max; end
        
        if ((propagAmplDeriv(hs, energy) * propagAmplDeriv(hs, energy + energy_step) < 0) 
            && propagAmplDeriv(hs, energy) < 0 )

            push!(hs.energy_levels, find_zero(zerofunc, energy, energy + energy_step)) 
        end

        if length(hs.energy_levels) == number_of_states
            break
        end

        energy      = energy_next
        energy_next = energy_next + energy_step
    end

    return hs.energy_levels
end


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function eigenenergy(hs::Heterostructure)
    eigenenergy(hs::Heterostructure, ERG_CALC_STEP) 
end
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function eigenenergy(hs::Heterostructure, step::Float64)
    eigenenergy(hs, step, 1000)
end
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function eigenenergy(hs::Heterostructure, number_of_states::Int=1000)
    eigenenergy(hs, ERG_CALC_STEP, number_of_states)
end
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -








# Calculate array of wavefunction coefficients
function waveFunctionCoefs(hs::Heterostructure, energy::Float64)

    coefs = Array(Array{Complex{Float64},1}, length(hs.layers))
    coefs[1] = [0.0; 1.0]       # Coefficients for the first layer A0=0, B0=1
    for i = 2 : length(hs.layers) 
        coefs[i] = propagationMatrix(hs, i-1, energy) * coefs[i-1]
    end    
    coefs[end][2] = complex(0)  # There must be zero, BN=0
    
    return coefs
end










# The Heterostructure x axis array or wavefunctions
function x_axis(hs::Heterostructure, density::Float64=DENSITY_OF_X) 

    lnum  = length(hs.layers)
    width = 0

    for i = 1:lnum
        width += hs.layers[i].width
    end
    
    density = (density > lnum / width) ? density : DENSITY_OF_X

    number = int(width * density)
    pfirst = hs.layers[ 1 ].center - 0.5*hs.layers[ 1 ].width
    plast  = hs.layers[end].center + 0.5*hs.layers[end].width

    # println("pfirst: $pfirst, plast: $plast, number: $number")
    # return [x for x in pfirst:step:plast]
    return linspace(pfirst, plast, number)
end








# Wave function values over passed x axis array
function waveFunction(  hs      ::Heterostructure, 
                        energy  ::Float64,
                        x       ::Array{Float64,1}  )

    len   = length(x)               # Array size
    y     = Array(Complex128, len)  # Y allocation
    coefs = waveFunctionCoefs(hs, energy)

    i     = 1                       # Coordinate index
    li    = 1                       # Layers index
    shift = hs.layers[1].center     # Coordinate shift
    wnum = waveNumber(hs.layers[1], energy)

    while i <= len
        if x[i] > shift + 0.5*hs.layers[li].width
            li += 1
            shift = hs.layers[li].center
            wnum = waveNumber(hs.layers[li], energy)
        end
        arg  = im * wnum * (x[i] - shift)
        y[i] = coefs[li][1] * exp(arg) + coefs[li][2] * exp(-arg)
        i += 1
    end

    area = numpy.trapz( abs(y.*y), x) 
    y /= sqrt(area) * WAVE_FUNCTION_DIVIDER

    return y
end






# Potential profile
function potentialForm(hs::Heterostructure)

    x = Array(Float64, 2*length(hs.layers))
    y = Array(Float64, 2*length(hs.layers))

    y_shift = -minimum(y)

    for i = 2 : 2 : 2*length(hs.layers)
        x[i-1] = hs.layers[int(i/2)].center - 0.5*hs.layers[int(i/2)].width
        x[ i ] = hs.layers[int(i/2)].center + 0.5*hs.layers[int(i/2)].width
        y[i-1] = hs.layers[int(i/2)].potential
        y[ i ] = hs.layers[int(i/2)].potential
    end
    return (x,y)
end






# TODO rename drawwf 
function drawwf(hs::Heterostructure, wfnum=nothing)
    # TODO Automatic wavefunction scaling depending on number of curves
    # to be shown and potential depth, avoid usage of WAVE_FUNCTION_DIVIDER
    # TODO Automatic scaling of the plot min/max of x and x
    # TODO Coloring and filling
    # TODO Title of the plot and axes
    # TODO Legend
    cm2nm = 1e7; # to scale in nm instead of 1e-7 cm

    if isempty(hs.energy_levels)
        number = (wfnum == nothing) ? 1000 : wfnum
        # Resolution 1 meV
        energy_levels = eigenenergy(hs, convert_ev2erg(0.001), number) 
    end

    if wfnum == nothing 
        wfnum = length(hs.energy_levels)
    end

    figure("Wave functions")

    band_x, band_y = potentialForm(hs)
    plot(band_x .* cm2nm, convert_erg2ev(band_y))
    
    wf_x = x_axis(hs)
    wf_x_nm = wf_x * cm2nm
    for i=1:wfnum
        wf_y = real(waveFunction(hs, hs.energy_levels[i], wf_x))
        plot( wf_x_nm, convert_erg2ev(hs.energy_levels[i] .+ wf_y) )
        plot( wf_x_nm[1:length(wf_x)-1:end], 
            convert_erg2ev([ hs.energy_levels[i]; hs.energy_levels[i]]) )
    end
end










# ==========================================================================
# TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST 

function test(roundto=3)

    # FIXME FIXME FIXME must be three states, one is missing!
    # FIXME           0.051483 eV
    # FIXME           0.198066 eV
    # FIXME           0.388992 eV -- missing
    # 
    # p = convert_ev2erg(0.4)     # potential, eV -> erg
    # ww = 8e-7                   # quantum well width, 8 nm
    # wb = 2e-5                   # barrier width, 20 nm
    # m = 0.067                   # particle mass



    p = convert_ev2erg(0.4)     # potential, eV -> erg
    # ww = 3e-7                   # quantum well width, 3 nm
    ww = 8e-7                   # quantum well width, 8 nm
    wb = 2e-5                   # barrier width, 20 nm
    m = 0.067                   # particle mass

    hs = Heterostructure([  Region(p, wb, -(ww+wb)/2, m) 
                            Region(0, ww,   0,        m)
                            Region(p, wb,  (ww+wb)/2, m) ], [])
    tic()
    energies = eigenenergy(hs, convert_ev2erg(0.001), 9) # Resolution 1 meV
    toc()

    energies_ev = convert_erg2ev(energies)

    println( join( ["$i: $(round(energies_ev[i],roundto)) eV" 
                    for i = 1:length(energies_ev)] ,"\n") )
    drawwf(hs)
end



end # of module tmm

