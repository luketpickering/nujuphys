import Base
import LinearAlgebra

abstract type PhysVector end

struct XYZVector <: PhysVector
    x::Real
    y::Real
    z::Real
end

struct XYZTVector <: PhysVector
    x::Real
    y::Real
    z::Real
    t::Real
end

abstract type PositionVector <: PhysVector end
struct Position3Vector <: PositionVector
    x::Real
    y::Real
    z::Real
end

struct Position4Vector <: PositionVector
    x::Real
    y::Real
    z::Real
    t::Real
end

abstract type DisplacementVector <: PhysVector end
struct Displacement3Vector <: DisplacementVector
    dx::Real
    dy::Real
    dz::Real
end

struct Displacement4Vector <: DisplacementVector
    dx::Real
    dy::Real
    dz::Real
    dt::Real
end

abstract type MomentumVector <: PhysVector end
struct Momentum3Vector <: MomentumVector
    px::Real
    py::Real
    pz::Real
end

struct Momentum4Vector <: MomentumVector
    px::Real
    py::Real
    pz::Real
    e::Real
end

ThreeVectors = Union{XYZVector,Position3Vector,Displacement3Vector,Momentum3Vector}
RotateableVectors = Union{XYZVector,Displacement3Vector,Momentum3Vector}
FourVectors = Union{XYZTVector,Position4Vector,Displacement4Vector,Momentum4Vector}

_threecomp(v::XYZTVector) = XYZVector(v.x,v.y,v.z)
_threecomp(v::Position4Vector) = Position3Vector(v.x,v.y,v.z)
_threecomp(v::Displacement4Vector) = Displacement3Vector(v.dx,v.dy,v.dz)
_threecomp(v::Momentum4Vector) = Momentum3Vector(v.px,v.py,v.pz)

_toarr(v::XYZVector) = [v.x,v.y,v.z]
_toarr(v::XYZTVector) = [v.x,v.y,v.z,v.t]

_toarr(v::Position3Vector) = [v.x,v.y,v.z]
_toarr(v::Position4Vector) = [v.x,v.y,v.z,v.t]

_toarr(v::Displacement3Vector) = [v.dx,v.dy,v.dz]
_toarr(v::Displacement4Vector) = [v.dx,v.dy,v.dz,v.dt]

_toarr(v::Momentum3Vector) = [v.px,v.py,v.pz]
_toarr(v::Momentum4Vector) = [v.px,v.py,v.pz,v.e]

_xyzfromarr(a) = XYZVector(a[1],a[2],a[3])
_xyztfromarr(a) = XYZVector(a[1],a[2],a[3],a[4])

mag2(v::ThreeVectors) = abs(_toarr(vectscale))
mag(v::ThreeVectors) = sqrt(mag2(v))

threemag2(v::FourVectors) = mag2(_threecomp(v))
threemag(v::FourVectors) = sqrt(threemag2(v))

mass2(v::Momentum4Vector) = (v.e^2 - (v.px^2 + v.py^2 + v.pz^2))
mass(v::Momentum4Vector) = sqrt(mass2(v))

scale(v::Position3Vector, xs::Real, ys::Real, zs::Real) = Position3Vector(v.x*xs,v.y*ys,v.z*zs)
scale(v::Displacement3Vector, xs::Real, ys::Real, zs::Real ) = Displacement3Vector(v.dx*xs,v.dy*ys,v.dz*zs)
scale(v::Momentum3Vector, xs::Real, ys::Real, zs::Real ) = Momentum3Vector(v.px*xs,v.py*ys,v.pz*zs)

vectscale(v::ThreeVectors, s::Real) = scale(v,s,s,s)

unit(v::RotateableVectors) = vectscale(v, 1.0/mag(v))

function rotxmat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  1  0  0  ;
       0  c -s  ;
       0  s  c  ]
end

function rotymat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  c  0  s  ;
       0  1  0  ;
      -s  0  c  ]
end

function rotzmat(ang) 
    c = cos(ang)
    s = sin(ang)
    [  c -s  0  ;
       s  c  0  ;
       0  0  1  ]
end

rotateaboutx(v::RotateableVectors, ang::Real) = rotxmat(ang) * _toarr(v)
rotateabouty(v::RotateableVectors, ang::Real) = rotymat(ang) * _toarr(v)
rotateaboutz(v::RotateableVectors, ang::Real) = rotzmat(ang) * _toarr(v)

function rotate(v::RotateableVectors, angyaw::Real, angpitch::Real, angroll::Real)
    rotzmat(angyaw) * rotymat(angpitch) * rotxmat(angroll) * _toarr(v)
end

function rotaboutaxismat(axis::RotateableVectors, ang::Real)
    # From https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    v = _toarr(unit(axis))

    c = cos(ang)
    s = sin(ang)
    diagc = LinearAlgebra.Diagonal([c,c,c])
    crossprodmat = [     0 -v[3]  v[2] ;
                      v[3]     0 -v[1] ;
                     -v[2]  v[1]     0 ]
    outerprod = v * v'

    diagc + s*crossprodmat + (1-c)*outerprod
end

function boost(r::T, dir::RotateableVectors, β::Real) where {T<:FourVectors}
    # From https://en.wikipedia.org/wiki/Lorentz_transformation#Proper_transformations
    n = _xyzfromarr(_toarr(unit(dir)))
    γ = 1.0/(sqrt(1.0 - β^2))

    boosted = 
    [ 1+(γ-1)*(n.x^2)    (γ-1)*n.x*n.y    (γ-1)*n.x*n.z  -γ*β*n.x;
        (γ-1)*n.y*n.z  1+(γ-1)*(n.y^2)    (γ-1)*n.y*n.z  -γ*β*n.y;
        (γ-1)*n.x*n.z    (γ-1)*n.y*n.z  1+(γ-1)*(n.z^2)  -γ*β*n.z;
             -γ*β*n.x         -γ*β*n.y         -γ*β*n.z         γ ] * _toarr(r)
    T(boosted[1],boosted[2],boosted[3],boosted[4])
end

function boost(r::T, v::RotateableVectors) where {T<:FourVectors}
    boost(r, v, mag(v))
end

rotateaboutaxis(v::RotateableVectors, axis::RotateableVectors, ang::Real) = rotaboutaxismat(axis,ang) * _toarr(v)

β(v::Momentum4Vector) = threemag(v)/v.e
γ(v::Momentum4Vector) = 1.0/sqrt(1.0 - β(v)^2)

function getboostvector(v::Momentum4Vector) 
  β(v) * unit(_threecomp(v))
end

LinearAlgebra.:⋅(v1::ThreeVectors, v2::ThreeVectors) = _toarr(v1) ⋅ _toarr(v2)
LinearAlgebra.:×(v1::ThreeVectors, v2::ThreeVectors) = _toarr(v1) × _toarr(v2)

#pos/pos = disp
Base.:-(v1::Position3Vector, v2::Position3Vector) = Displacement3Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)
Base.:-(v1::Position4Vector, v2::Position4Vector) = Displacement4Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.t - v2.t)

#disp/disp = disp
Base.:+(v1::Displacement3Vector, v2::Displacement3Vector) = Displacement3Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz)
Base.:+(v1::Displacement4Vector, v2::Displacement4Vector) = Displacement4Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz, v1.dt + v2.dt)
Base.:-(v1::Displacement3Vector, v2::Displacement3Vector) = Displacement3Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz)
Base.:-(v1::Displacement4Vector, v2::Displacement4Vector) = Displacement4Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz, v1.dt - v2.dt)

#mom/mom = mom
Base.:+(v1::Momentum3Vector, v2::Momentum3Vector) = Momentum3Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz)
Base.:+(v1::Momentum4Vector, v2::Momentum4Vector) = Momentum4Vector(v1.dx + v2.dx, v1.dy + v2.dy, v1.dz + v2.dz, v1.dt + v2.dt)
Base.:-(v1::Momentum3Vector, v2::Momentum3Vector) = Momentum3Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz)
Base.:-(v1::Momentum4Vector, v2::Momentum4Vector) = Momentum4Vector(v1.dx - v2.dx, v1.dy - v2.dy, v1.dz - v2.dz, v1.dt - v2.dt)

#pos/disp = pos
Base.:+(v1::Position3Vector, v2::Displacement3Vector) = Position3Vector(v1.x + v2.dx, v1.y + v2.dy, v1.z + v2.dz)
Base.:+(v1::Position4Vector, v2::Displacement4Vector) = Position4Vector(v1.x + v2.dx, v1.y + v2.dy, v1.z + v2.dz, v1.t + v2.dt)
Base.:-(v1::Position3Vector, v2::Displacement3Vector) = Position3Vector(v1.x - v2.dx, v1.y - v2.dy, v1.z - v2.dz)
Base.:-(v1::Position4Vector, v2::Displacement4Vector) = Position4Vector(v1.x - v2.dx, v1.y - v2.dy, v1.z - v2.dz, v1.t - v2.dt)

Base.:*(s::Real, v::ThreeVectors) = vectscale(v, s)
Base.:*(v::ThreeVectors, s::Real) = vectscale(v, s)
Base.:/(v::ThreeVectors, s::Real) = vectscale(v, 1.0/s)
Base.:-(v::ThreeVectors) = vectscale(v, -1)

Base.show(io::IO, v::XYZVector) = print(io, "[x: ",v.x,", y:",v.y,", z:",v.z,"]")
Base.show(io::IO, v::XYZTVector) = print(io, "[x: ",v.x,", y: ",v.y,", z:",v.z,", t:",v.t,"]")
Base.show(io::IO, v::Position3Vector) = print(io, "[pos: ",v.x,",",v.y,",",v.z,"]")
Base.show(io::IO, v::Position4Vector) = print(io, "[pos: ",v.x,",",v.y,",",v.z,",",v.t,"]")
Base.show(io::IO, v::Displacement3Vector) = print(io, "[disp: ",v.dx,",",v.dy,",",v.dz,"]")
Base.show(io::IO, v::Displacement4Vector) = print(io, "[disp: ",v.dx,",",v.dy,",",v.dz,",",v.dt,"]")
Base.show(io::IO, v::Momentum3Vector) = print(io, "[mom: ",v.px,",",v.py,",",v.pz,"]")
Base.show(io::IO, v::Momentum4Vector) = print(io, "[mom: ",v.px,",",v.py,",",v.pz,",",v.e,"]")
