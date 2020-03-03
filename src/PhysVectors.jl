module PhysVectors

import Base: convert
import LinearAlgebra

abstract type PhysVector end
abstract type PositionVector <: PhysVector end
abstract type DisplacementVector <: PhysVector end
abstract type MomentumVector <: PhysVector end

"""
  Generic three vector, used where you don't want the type system to restrict operations.
"""
struct XYZVector <: PhysVector
    x::Real
    y::Real
    z::Real
end
Base.show(io::IO, v::XYZVector) = print(io, "[x: ",v.x,", y:",v.y,", z:",v.z,"]")
Base.:+(v1::XYZVector, v2::XYZVector) = XYZVector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)
Base.:-(v1::XYZVector, v2::XYZVector) = XYZVector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)

"""
  Generic four vector, used where you don't want the type system to restrict operations.
"""
struct XYZTVector <: PhysVector
    x::Real
    y::Real
    z::Real
    t::Real
end
Base.show(io::IO, v::XYZTVector) = print(io, "[x: ",v.x,", y: ",v.y,", z:",v.z,", t:",v.t,"]")
Base.:+(v1::XYZTVector, v2::XYZTVector) = XYZTVector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.t + v2.t)
Base.:-(v1::XYZTVector, v2::XYZTVector) = XYZTVector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.t - v2.t)

"""
    Three vector describing a spatial displacement.

    Displacements can be summed, subtracted, scaled, and rotated.
"""
struct Displacement3Vector <: DisplacementVector
    dx::Real
    dy::Real
    dz::Real
end
Base.show(io::IO, v::Displacement3Vector) = print(io, "[disp: ",v.dx,",",v.dy,",",v.dz,"]")
Base.:+(v1::Displacement3Vector, v2::Displacement3Vector) = Displacement3Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz)
Base.:-(v1::Displacement3Vector, v2::Displacement3Vector) = Displacement3Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz)

"""
    Four vector describing a spacetime displacement.

    Four displacements can be boosted, summed, subtracted, and their spatial components can be scaled and rotated.
"""
struct Displacement4Vector <: DisplacementVector
    dx::Real
    dy::Real
    dz::Real
    dt::Real
end
Base.show(io::IO, v::Displacement4Vector) = print(io, "[disp: ",v.dx,",",v.dy,",",v.dz,",",v.dt,"]")
Base.:+(v1::Displacement4Vector, v2::Displacement4Vector) = Displacement4Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz, v1.dt + v2.dt)
Base.:-(v1::Displacement4Vector, v2::Displacement4Vector) = Displacement4Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz, v1.dt - v2.dt)

"""
    Three vector describing a spatial position.

    Can be handled by a restricted set of operators, e.g. you cannot add positions, but you can subtract them, producing a displacement vector.

    Positions can be rotated and translated (with a displacement vector).
"""
struct Position3Vector <: PositionVector
    x::Real
    y::Real
    z::Real
end
Base.show(io::IO, v::Position3Vector) = print(io, "[pos: ",v.x,",",v.y,",",v.z,"]")

Base.:-(v1::Position3Vector, v2::Position3Vector) = Displacement3Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)
Base.:+(v1::Position3Vector, v2::Displacement3Vector) = Position3Vector(v1.x + v2.dx, v1.y + v2.dy, v1.z + v2.dz)
Base.:-(v1::Position3Vector, v2::Displacement3Vector) = Position3Vector(v1.x - v2.dx, v1.y - v2.dy, v1.z - v2.dz)

"""
    Four vector describing a spacetime position.

    Can be handled by a restricted set of operators, e.g. you cannot add positions, but you can subtract them, producing a displacement vector.

    Four positions can be boosted, and translated, their spatial components can be rotated.
"""
struct Position4Vector <: PositionVector
    x::Real
    y::Real
    z::Real
    t::Real
end
Base.show(io::IO, v::Position4Vector) = print(io, "[pos: ",v.x,",",v.y,",",v.z,",",v.t,"]")

Base.:-(v1::Position4Vector, v2::Position4Vector) = Displacement4Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.t - v2.t)
Base.:+(v1::Position4Vector, v2::Displacement4Vector) = Position4Vector(v1.x + v2.dx, v1.y + v2.dy, v1.z + v2.dz, v1.t + v2.dt)
Base.:-(v1::Position4Vector, v2::Displacement4Vector) = Position4Vector(v1.x - v2.dx, v1.y - v2.dy, v1.z - v2.dz, v1.t - v2.dt)

"""
    Three vector describing a momentum.

    Three momenta can be summed, subtracted, scaled, and rotated.
"""
struct Momentum3Vector <: MomentumVector
    px::Real
    py::Real
    pz::Real
end
Base.show(io::IO, v::Momentum3Vector) = print(io, "[mom: ",v.px,",",v.py,",",v.pz,"]")
Base.:+(v1::Momentum3Vector, v2::Momentum3Vector) = Momentum3Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz)
Base.:-(v1::Momentum3Vector, v2::Momentum3Vector) = Momentum3Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz)

"""
    Four vector describing a four momentum.

    Four momenta can be boosted, summed, subtracted.
"""
struct Momentum4Vector <: MomentumVector
    px::Real
    py::Real
    pz::Real
    e::Real
end
Base.show(io::IO, v::Momentum4Vector) = print(io, "[mom: ",v.px,",",v.py,",",v.pz,",",v.e,"]")
Base.:+(v1::Momentum4Vector, v2::Momentum4Vector) = Momentum4Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz, v1.dt + v2.dt)
Base.:-(v1::Momentum4Vector, v2::Momentum4Vector) = Momentum4Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz, v1.dt - v2.dt)

"""
    Internal function for getting the spatial component of a fourvector while respecting the vector type hierarchy
"""
function _threecomp(v::XYZTVector) 
    XYZVector(v.x,v.y,v.z)
end
_threecomp(v::Position4Vector) = Position3Vector(v.x,v.y,v.z)
_threecomp(v::Displacement4Vector) = Displacement3Vector(v.dx,v.dy,v.dz)
_threecomp(v::Momentum4Vector) = Momentum3Vector(v.px,v.py,v.pz)

"""
    Internal function for array-izing typed vectors to use julia BLAS
"""
function _toarr(v::XYZVector)
    [v.x,v.y,v.z]
end
_toarr(v::XYZTVector) = [v.x,v.y,v.z,v.t]
_toarr(v::Position3Vector) = [v.x,v.y,v.z]
_toarr(v::Position4Vector) = [v.x,v.y,v.z,v.t]
_toarr(v::Displacement3Vector) = [v.dx,v.dy,v.dz]
_toarr(v::Displacement4Vector) = [v.dx,v.dy,v.dz,v.dt]
_toarr(v::Momentum3Vector) = [v.px,v.py,v.pz]
_toarr(v::Momentum4Vector) = [v.px,v.py,v.pz,v.e]

ThreeVectors = Union{ 
    XYZVector,
    Position3Vector,
    Displacement3Vector,
    Momentum3Vector
}

LinearAlgebra.:⋅(v1::ThreeVectors, v2::ThreeVectors) = _toarr(v1) ⋅ _toarr(v2)
LinearAlgebra.:×(v1::ThreeVectors, v2::ThreeVectors) = _toarr(v1) × _toarr(v2)

mag2(v::ThreeVectors) = LinearAlgebra.norm(_toarr(v))
mag(v::ThreeVectors) = sqrt(mag2(v))

DirectionVectors = Union{ 
    XYZVector,
    Displacement3Vector,
    Momentum3Vector
}

Base.convert(::Type{XYZVector}, x::Displacement3Vector) = XYZVector(x.dx, x.dy, x.dz)
Base.convert(::Type{XYZVector}, x::Momentum3Vector) = XYZVector(x.px, x.py, x.pz)
Base.convert(::Type{XYZTVector}, x::Displacement4Vector) = XYZTVector(x.dx, x.dy, x.dz, x.dt)
Base.convert(::Type{XYZTVector}, x::Momentum4Vector) = XYZTVector(x.px, x.py, x.pz, x.e)


"""
    Build a three rotation matrix around the x axis
"""
function rotxmat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  1  0  0  ;
       0  c -s  ;
       0  s  c  ]
end

"""
    Build a three rotation matrix around the y axis
"""
function rotymat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  c  0  s  ;
       0  1  0  ;
      -s  0  c  ]
end

"""
    Build a three rotation matrix around the z axis
"""
function rotzmat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  c -s  0  ;
       s  c  0  ;
       0  0  1  ]
end

"""
    Rotate a three vector by rotation matrix, respects the physics vector type
"""
function rotate(v::T, matrix::Array{R,2}) where {T<:ThreeVectors, R<:Real}
    rotated = matrix * _toarr(v)
    T(rotated[1], rotated[2], rotated[3])
end

"""
    Rotate a three vector around x by an angle, respects the physics vector type
"""
rotateaboutx(v::T, ang::Real) where {T<:ThreeVectors} = rotate(v, rotxmat(ang))
"""
    Rotate a three vector around y by an angle, respects the physics vector type
"""
rotateabouty(v::T, ang::Real) where {T<:ThreeVectors} = rotate(v, rotymat(ang))
"""
    Rotate a three vector around z by an angle, respects the physics vector type
"""
rotateaboutz(v::T, ang::Real) where {T<:ThreeVectors} = rotate(v, rotzmat(ang))

"""
    Rotate a three vector by roll around x, then by pitch around y, then by yaw around z, respects the physics vector type
"""
function rotatexyz(v::T, angyaw::Real, angpitch::Real, angroll::Real)  where {T<:ThreeVectors}
    rotate(v, rotzmat(angyaw) * rotymat(angpitch) * rotxmat(angroll))
end

"""
    Build a three rotation matrix around an arbitrary three axis
"""
function rotationmat(axis::DirectionVectors, ang::Real)
    # From https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    v = _toarr(unit(axis))

    c = cos(ang)
    s = sin(ang)
    diagc = LinearAlgebra.Diagonal([c,c,c])
    crossprodmat = [     0 -v[3]  v[2] ;
                      v[3]     0 -v[1] ;
                     -v[2]  v[1]     0 ]
    outerprod = v * v'

    diagc + s * crossprodmat + (1-c) * outerprod
end

"""
    Rotate a three vector by roll around x, then by pitch around y, then by yaw around z, respects the physics vector type
"""
rotateaboutaxis(v::T, axis::DirectionVectors, ang::Real) where {T<:ThreeVectors} = rotate(v, rotationmat(axis, ang))

ScaleableVectors = Union{ 
    XYZVector,
    Displacement3Vector,
    Momentum3Vector
}

"""
    Scales a three vector by a different scale along each principal axis.
"""
function scale(v::T, xs::Real, ys::Real, zs::Real) where {T<:ScaleableVectors} 
    scaled = [xs, ys, zs] .* _toarr(v)
    T(scaled[1], scaled[2], scaled[3])
end

"""
    Scales a three vector by a single scalefactor
"""
scale(v::ScaleableVectors, s::Real) = scale(v,s,s,s)

"""
    Normalizes a vector to length 1
"""
unit(v::ScaleableVectors) = scale(v, 1.0/mag(v))

function Base.:*(s::Real, v::T) where {T<:ScaleableVectors} 
    scale(v, s)
end
Base.:*(v::T, s::Real) where {T<:ScaleableVectors} = s * v
Base.:/(v::T, s::Real) where {T<:ScaleableVectors} = (1.0/s) * v
Base.:-(v::T) where {T<:ScaleableVectors} = -1.0 * v

FourVectors = Union{ 
    XYZTVector,
    Position4Vector,
    Displacement4Vector,
    Momentum4Vector
}

function LinearAlgebra.:⋅(v1::FourVectors, v2::FourVectors)
    a1 = _toarr(v1)
    a2 = _toarr(v2)
    (a1[4] * a2[4]) - LinearAlgebra.:⋅(a1[1:3], a2[1:3])
end

threemag2(v::FourVectors) = mag2(_threecomp(v))
threemag(v::FourVectors) = sqrt(threemag2(v))

function fourmag2(v::FourVectors)
    a = _toarr(v)
    a[4]^2 - (a[1]^2 + a[2]^2 + a[3]^2)
end

function fourmag(v::FourVectors)
    sqrt(fourmag2(v))
end

BoostableVectors = Union{
    XYZTVector,
    Position4Vector,
    Momentum4Vector
}

function boost(r::T, dir::DirectionVectors, β::Real) where {T<:BoostableVectors}
    # From https://en.wikipedia.org/wiki/Lorentz_transformation#Proper_transformations
    n = Base.convert(XYZVector, unit(dir))
    γ = 1.0/(sqrt(1.0 - β^2))

    boosted = 
    [ 1+(γ-1)*(n.x^2)    (γ-1)*n.x*n.y    (γ-1)*n.x*n.z  -γ*β*n.x;
        (γ-1)*n.y*n.z  1+(γ-1)*(n.y^2)    (γ-1)*n.y*n.z  -γ*β*n.y;
        (γ-1)*n.x*n.z    (γ-1)*n.y*n.z  1+(γ-1)*(n.z^2)  -γ*β*n.z;
             -γ*β*n.x         -γ*β*n.y         -γ*β*n.z         γ ] * _toarr(r)
    T(boosted[1],boosted[2],boosted[3],boosted[4])
end

function boost(r::T, v::XYZVector) where {T<:BoostableVectors}
    boost(r, v, mag(v))
end

β(v::Momentum4Vector) = threemag(v)/v.e
γ(v::Momentum4Vector) = 1.0/sqrt(1.0 - β(v)^2)

function getboostvector(v::Momentum4Vector) 
  β(v) * unit(_threecomp(v))
end

mass2(v::Momentum4Vector) = fourmag2(v)
mass(v::Momentum4Vector) = fourmag(v)

end