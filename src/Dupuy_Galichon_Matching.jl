module Dupuy_Galichon_Matching 

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using LinearAlgebra


"""
    convertion(x::Any) 

Converts element to Float64, usefull to convert String. It is especially coded to convert the strings using commas. \\
convertion("0,01") --> returns Float64(0.01)
"""
function convertion(x::Any) 
    if x isa String
        if occursin(",", x) == true
             c = 1
             xx = x[c]
                 while xx != ','
                     c += 1
                    xx = x[c]
                  end
        x = string(x[1:(c-1)],".",x[(c+1):length(x)])
        x = parse(Float64, x)
        else
        x = parse(Float64, x)
        end
    else
    convert(Float64, x) 
    end
end


"""
    convertion_df(x::DataFrame) 

Converts every elements of a DataFrame to type Float64. \\
Based on convertion()
"""
function convertion_df!(x::DataFrame)
    for i in 1:size(x)[2]
        x[i] = convertion.(x[i])
    end      
end

"""
    z_diag_f(x::Matrix, f::Function)

Mimick the matlab function diag(). \\
Apply a function on a Matrix columns (such as mean() or var()). \\
Returns a zero matrix with the computed values in its diagonal.
"""
function z_diag_f(x::Matrix, f::Function)
    z = zeros(size(x)[2],size(x)[2]) 
    f_x = f.(eachcol(x))
    z[diagind(z)] = f_x
    return(z)
end

"""
    z_diag(x::Vector)

Returns a zero matrix with the input vector as diagonal.
"""
function z_diag(x::Vector)
    z = zeros(length(x),length(x))
    z[diagind(z)] = x
    return(z)
end

end