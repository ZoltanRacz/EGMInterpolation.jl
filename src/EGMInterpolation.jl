module EGMInterpolation

using DocStringExtensions: FIELDS, TYPEDSIGNATURES, TYPEDEF # For easier documentation. 

export EGMInterpolatedFunction,
       evaluate

include("interpolation.jl")

end
