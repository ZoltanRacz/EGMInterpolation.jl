"""
gives the index of the higher z point which will be taken into account when interpolating
"""
function findz(zgrid::AbstractVector,z::Real)
    zdim = length(zgrid)
    if z<zgrid[1]
    # if z<=zgrid[1]
        zhin = 2
    elseif z>=zgrid[end]
    # elseif z>zgrid[end]
        zhin = zdim
    else
        zhin::Int64 = findfirst(x->x>z,zgrid)
        # zhin = searchsortedfirst(zgrid, z)
    end
    return zhin
end

"""
gives the index of the higher a point which will be taken into account when interpolating, given the relevant z index
"""
function finda(agrid::AbstractVector,zgrid::AbstractVector,zin::Integer,a::Real)
    #amin = tg.amins[zin]
    zdim = length(zgrid)
    adim = length(agrid[zin])
    if a<agrid[zin][1]
    # if a<=agrid[zin][1]
        ahin = 2
    elseif a>=agrid[zin][adim]
    # elseif a>agrid[zin][adim]
        ahin = adim
    else
#        ahin = findfirst(x->x>a,@view agrid[zin])::Int64
        ahin::Int64 = findfirst(x->x>a,agrid[zin])
        # ahin = searchsortedfirst(agrid[zin], a)
    end
    return ahin
end

"""

linear inter/extrapolation from the values at the 4 nearby (z,a) grid points
"""
function lin_trap(a,z,al1,ah1,al2,ah2,zl,zh,vl1,vh1,vl2,vh2)
    zrat = (z-zl)/(zh-zl)
    arat1 = (a-al1)/(ah1-al1)
    arat2 = (a-al2)/(ah2-al2)
    intval = (1-zrat)*(arat1*vh1+(1-arat1)*vl1) + zrat*(arat2*vh2+(1-arat2)*vl2)
    return intval
end

"""

Perform 2D linear interpolation at (z,a) given grids and values, with arbitrary z.

agrid and vals are both vectors of vectors. This is necessary as coh grids depends on permannet income state. For example, agrid[zi] gives the grid for coh corresponding to zgrid[zi]. 
"""
function evaluate(agrid::AbstractVector,vals::AbstractVector,zgrid::AbstractVector,z::Real,a::Real)
    zhin = findz(zgrid,z)
    zlin = zhin-1
    ah1in = finda(agrid,zgrid,zlin,a)
    ah2in = finda(agrid,zgrid,zhin,a)
    al1in = ah1in-1
    al2in = ah2in-1
    zl = zgrid[zlin]
    zh = zgrid[zhin]
    al1 = agrid[zlin][al1in]
    ah1 = agrid[zlin][ah1in]
    al2 = agrid[zhin][al2in]
    ah2 = agrid[zhin][ah2in]
    return lin_trap(a,z,al1,ah1,al2,ah2,zl,zh,
    vals[zlin][al1in],vals[zlin][ah1in],vals[zhin][al2in],vals[zhin][ah2in])
end

"""

Performs 2D linear interpolation at (z,a) given grids and values, with z on the grid, given with index (integer).

agrid and vals are both vectors of vectors. This is necessary as coh grids depends on permannet income state. For example, agrid[zi] gives the grid for coh corresponding to zgrid[zi]. 
"""
function evaluate(agrid::AbstractVector,vals::AbstractVector,zgrid::AbstractVector,zi::Integer,a::Real)

    ahin = finda(agrid,zgrid,zi,a)
    alin = ahin-1
    return vals[zi][ahin] - (a-agrid[zi][ahin])/(agrid[zi][alin]-agrid[zi][ahin])*(vals[zi][ahin]-vals[zi][alin])
end

#=
From here: unused version. Might be valuable, but for some reason we chose the simpler method. Comment out or delete, but save the code somewhere!
=#

function finda2(agrid::AbstractVector,zgrid::AbstractVector,zin::Integer,a::Real)
    @warn("not updated")
    #amin = tg.amins[zin]
    zdim = length(zgrid)
    adim = length(agrid[zin])
    if a>=agrid[zin][adim]
        ahin = adim
    else
#        ahin = findfirst(x->x>a,@view agrid[zin])::Int64
        ahin = findfirst(x->x>a,agrid[zin])::Int64
    end
    return ahin
end

"""
linear inter/extrapolation from 3 values nearby (z,a) grid points. could be used
for more refined version taking a better care of kinks from the boundary constraint
and the fixed participation cost

currently unused, who knows why
"""
function lin_tri(a,z,asolo,al,ah,zsolo,zd,vsolo,vl,vh)
    if z==zsolo
        return vsolo
    else
        zrat = (z-zsolo)/(zd-zsolo)
        arat = (a-(zrat*al+(1-zrat)*asolo))/(zrat*(ah-al))
        intval = (1-zrat)*vsolo + zrat*(arat*vh+(1-arat)*vl)
        return intval
    end
end

"""
version of evaluate based on `lin_tri`
"""
function evaluate2(agrid::AbstractVector,vals::AbstractVector,zgrid::AbstractVector,z::Real,a::Real)
    zhin = findz(tg,z)
    zlin = zhin-1
    ah1in = finda2(agrid,tg,zlin,a)
    ah2in = finda2(agrid,tg,zhin,a)
    #ah1in == 1 && ah2in == 1 && @error("a is smaller than both minima")
    if ah1in == 1 && ah2in == 1
        return 0.0
    end
    if ah2in==1
        return lin_tri(a,z,agrid[zhin][1],agrid[zlin][ah1in-1],agrid[zlin][ah1in],
        zgrid[zhin],zgrid[zlin],vals[zhin][1],vals[zlin][ah1in-1],vals[zlin][ah1in])
    end
    if ah1in==1
        return lin_tri(a,z,agrid[zlin][1],agrid[zhin][ah2in-1],agrid[zhin][ah2in],
        zgrid[zlin],zgrid[zhin],vals[zlin][1],vals[zhin][ah2in-1],vals[zhin][ah2in])
    end
    al1in = ah1in-1
    al2in = ah2in-1
    zl = zgrid[zlin]
    zh = zgrid[zhin]
    al1 = agrid[zlin][al1in]
    ah1 = agrid[zlin][ah1in]
    al2 = agrid[zhin][al2in]
    ah2 = agrid[zhin][ah2in]
    zrat = (z-zl)/(zh-zl)
    if agrid[zlin][1]<tg.amins[zlin]
        b1 = 1
    else
        b1 = 2
    end
    if agrid[zhin][1]<tg.amins[zhin]
        b2 = 1
    else
        b2 = 2
    end
    if ah1in>b1 && ah2in>b2
        lin_trap(a,z,al1,ah1,al2,ah2,zl,zh,
        vals[zlin][al1in],vals[zlin][ah1in],vals[zhin][al2in],vals[zhin][ah2in])
    elseif ah1in<=b1 && ah2in<=b2
        lin_trap(a,z,al1,ah1,al2,ah2,zl,zh,
        vals[zlin][al1in],vals[zlin][ah1in],vals[zhin][al2in],vals[zhin][ah2in])
    elseif ah1in<=b1 && ah2in>b2
        if a > zrat*agrid[zhin][b2] + (1-zrat)*ah1
            return lin_tri(a,z,ah1,al2,ah2,zl,zh,vals[zlin][ah1in],vals[zhin][al2in],vals[zhin][ah2in])
        else
            return lin_tri(a,z,al2,al1,ah1,zh,zl,vals[zhin][al2in],vals[zlin][al1in],vals[zlin][ahlin])
        end
    elseif ah1in>b1 && ah2in<=b2
        if a > zrat*ah2 + (1-zrat)*agrid[zlin][b1]
            return lin_tri(a,z,ah2,al1,ah1,zh,zl,vals[zhin][ah2in],vals[zlin][al1in],vals[zlin][ah1in])
        else
            return lin_tri(a,z,al1,al2,ah2,zl,zh,vals[zlin][al1in],vals[zhin][al2in],vals[zhin][ah2in])
        end
    else
        @error("what then?")
    end
end

"""
$(TYPEDEF)
# Description
Interpolation object
# Fields
$(FIELDS)
"""
struct EGMInterpolatedFunction{T<:AbstractFloat}
    "Shared grid, for uniform dimension"
    zgrid::Vector{T}
    "Different grids for each z point"
    agrids::Vector{Vector{T}}
    "Function values at grid points"
    vals::Vector{Vector{T}}
end

evaluate(egmif::EGMInterpolatedFunction,z::Real,a::Real) = evaluate(egmif.agrids,egmif.vals,egmif.zgrid,z,a)

(egmif::EGMInterpolatedFunction)(z::Real,a::Real) = evaluate(egmif,z,a)