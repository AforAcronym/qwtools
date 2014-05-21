
export varshni
function varshni(e0::Float64, a::Float64, b::Float64, t::Float64) 
    return e0 - a * t * t / (t + b)
end

export vegard
function vegard(x::Float64, p1::Float64, p2::Float64, bow::Float64) 
    return p1 * x + p2 * (1 - x) - bow * x * (1 - x)
end




export AlGaN_gap
function AlGaN_gap(t::Float64, x::Float64)
    E_GaN = 3.42 # at T = 0
    E_AlN = 6.08 # at T = 0
    e_bow = 1.1
    
    a_AlN = 2.63e-3
    a_GaN = 0.94e-3
    a_bow = 2.15e-3
    
    b_AlN = 2082.0 
    b_GaN =  791.0
    b_bow = 1561.0
    
    egap0 = vegard(x, E_AlN, E_GaN, e_bow)
    alpha = vegard(x, a_AlN, a_GaN, a_bow)
    beta  = vegard(x, b_AlN, b_GaN, b_bow)
    
    return varshni(egap0, alpha, beta, t)
end



export AlGaN_gap_tbl
function AlGaN_gap_tbl(t_step::Float64, x_step::Float64, t_max::Float64=300.0, print::Bool=false)
    if ! (0 < x_step <= 1)
        error("x_step must be in range (0,1]")
    end
    if t_step > t_max
        error("Temperature step cannot be greater than maximum: $(t_max)")
    end
    if t_step < 0
        error("Temperature step cannot be negative")
    end
    if t_max < 0
        error("Temperature maximum cannot be negative")
    end
    t_range = 0 : t_step : t_max
    x_range = 0 : x_step : 1

    tbl = Array(Float64, length(t_range), length(x_range))
    counter = 1
    for x = 0:x_step:1, t = 0:t_step:t_max
        tbl[counter] = AlGaN_gap(t, x)
        counter += 1
    end
    
    return tbl
end



export print_tbl_x
function print_tbl_x(tbl::Array{Float64,2}, acc::Int64=6)
    rows, cols = size(tbl)
    for i = 1:cols
        print_joined(STDOUT, tbl[:,i], "\t")
        println()
    end
end



export print_tbl_t
function print_tbl_t(tbl::Array{Float64,2}, acc::Int64=6)
    rows, cols = size(tbl)
    for i = 1:rows
        print_joined(STDOUT, tbl[i,:], "\t")
        println()
    end
end



export print_AlGaN_gap_tbl
function print_AlGaN_gap_tbl(t_step::Float64, x_step::Float64, t_max::Float64=300.0)
    tbl = AlGaN_gap_tbl(t_step, x_step, t_max)
    xrange = 0 : x_step : 1
    trange = 0 : t_step : t_max
    
    # for x = xrange
    #     print_joined(STDOUT, [x, [AlGaN_gap(t,x) for t in trange]], "\t")
    #     println()
    # end
    
    # println()
    # println()
    # println()
    
    print_joined(STDOUT, ["t \\ x", [x for x in xrange]], "\t")
    println()
    for t = trange
        print_joined(STDOUT, [t, [AlGaN_gap(t,x) for x in xrange]], "\t")
        println()
    end

end



# function S = algan(x)
# % Parameters set of Al(x)Ga(1-x)N alloy, 0<=x<=1

# xx  = [ 1-x x ]';
# xxx = [ xx' prod(xx) ]';


# % Alloy composition
# % S.x = x;


# % Energy gap, eV
# % 
# %        GaN     AlN     Bowing factor
# S.Eg = [ 3.42    6.08       1.1         ] * xxx;


# % Lattice constants, A
# %       GaN     AlN
# S.a = [ 3.18    3.11] * xx;
# S.c = [ 5.18    4.98] * xx;


# % Dielectric constant
# %          GaN   AlN
# S.diel = [ 8.9   8.5 ] * xx;


# % Piezoelectric constants, C/m^2
# %           GaN      AlN
# S.pz13 = [ -0.33    -0.58  ] * xx;
# S.pz33 = [  0.65     1.55  ] * xx;


# % Spontaneous polarization, C/m^2
# % 
# % Source:
# % Yu, Dang, Asbeck, Lau, Sullivan
# % J. Vac. Sci. Technol. B 17(4) 1742 Jul/Aug 1999
# % 
# %          GaN     AlN
# S.Psp = [ -0.029  -0.081  ] * xx;



# % Elastic constants, GPa
# % Data from www.ioffe.ru
# % 
# % Source:
# % http://www.ioffe.ru/SVA/NSM/Semicond/GaN/mechanic.html
# % http://www.ioffe.ru/SVA/NSM/Semicond/AlN/mechanic.html
# % 
# %         GaN    AlN
# S.c11 = [ 390    410  ] * xx;
# S.c12 = [ 145    149  ] * xx;
# S.c13 = [ 106     99  ] * xx;
# S.c33 = [ 398    389  ] * xx;
# % Not needed:
# % S.c44 = [ 105    125  ] * xx;


# S.G = 2 * ( S.c11 + S.c12 - 2*( S.c13^2 / S.c33 ));
# % S.D = 2 * ( S.c12 + 2*( S.c13^2 / S.c33 )) / ( S.c11 + 2*S.c12 );



# % Carriers masses, in units of free electron mass
# % m0 = 9.1095e-28 g
# % 
# % Source: 
# % http://www.ioffe.ru/SVA/NSM/Semicond/GaN/basic.html
# % http://www.ioffe.ru/SVA/NSM/Semicond/AlN/basic.html
# % 
# %           GaN    AlN
# S.me   = [  0.2    0.4   ] * xx;
# S.mhhx = [  1.6   10.42  ] * xx;
# S.mhhz = [  1.1    3.53  ] * xx;
# S.mlhx = [  1.1    0.24  ] * xx;
# S.mlhz = [  0.15   3.53  ] * xx;

