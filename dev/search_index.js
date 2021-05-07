var documenterSearchIndex = {"docs":
[{"location":"#Our-Replication-of-Dupuy-and-Galichon-\"Personality-Traits-and-the-Marriage-Market\"-(2014)","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"","category":"section"},{"location":"","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"@Gabriele Buontempo @mattjs07","category":"page"},{"location":"","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"This replication study was part of our evaluation for the course Numerical Methods at SciencesPo Paris in Spring 2021","category":"page"},{"location":"","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"In this replication study, we translate the original matlab code into Julia lang. You will find in this documentation the main function \" xxxxx \" as well as the intermediary functions we built for this function.","category":"page"},{"location":"","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"Modules = [Dupuy_Galichon_Matching]","category":"page"},{"location":"#Dupuy_Galichon_Matching.CovVarEdsi-Tuple{Matrix{T} where T, Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.CovVarEdsi","text":"CovVarEdsi(X::Matrix,Y::Matrix )\n\nReturns CovXYXY , CovXXXX, CovYYYY, CovXXYY\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.EmpCov-Tuple{Matrix{T} where T, Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.EmpCov","text":"EmpCov(X::Matrix,Y::Matrix)\n\nReturns the empirical covariance of two matrices.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.ObjectiveFunction-Tuple{Matrix{T} where T, Matrix{T} where T, Matrix{T} where T, Int64, Int64, Vararg{Any, N} where N}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.ObjectiveFunction","text":"ObjectiveFunction(X::Matrix,Y::Matrix,sigmaHat::Matrix,T::Int64,A... )\n\nObjective function used to compute the optimal Affinity Matrix. \nA is splat as it is the argument over which we arg min using JuMP. \ncolX is required to reshape A into a matrix.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.OptimalPi-Tuple{Matrix{T} where T, Matrix{T} where T, Matrix{T} where T, Int64}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.OptimalPi","text":"OptimalPi(auxVar::Matrix, X::Matrix, Y::Matrix, T::Int64)\n\nReturns the Optimal Pi, used in ObjectiveFunction() for optimization.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.StartingValues-Tuple{Matrix{T} where T, Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.StartingValues","text":"StartingValues(X::Matrix, Y::Matrix)  \n\nN,M = size(X)\nu0 = log(N) * ones(N,1)\nv0 = -ones(N,1)\nreturn u0, v0\n\nUsed in OptimalPi() function.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.convertion-Tuple{Any}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.convertion","text":"convertion(x::Any)\n\nConverts element to Float64, usefull to convert String. It is especially coded to convert the strings using commas. \nconvertion(\"0,01\") –> returns Float64(0.01)\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.convertion_df!-Tuple{DataFrames.DataFrame}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.convertion_df!","text":"convertion_df(x::DataFrame)\n\nConverts every elements of a DataFrame to type Float64. \nBased on convertion()\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.discretePdf-Tuple{Matrix{T} where T, Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.discretePdf","text":"discretePdf(X::Matrix, Y::Matrix)\n\nReturns two vectors of length equal to the matrices highest dimension (N). \nEach element equals to 1 divided by the number of rows. Used in OptimalPi() function.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.length_-Tuple{Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.length_","text":"length_(X::Matrix)\n\nMimick matlab length() function by returning the length of the biggest dimension.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.mm-Tuple{Matrix{T} where T, Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.mm","text":"mm(x::Matrix, y::Matrix)\n\nPerforms matrice multiplication using loops. \nUsefull for optimisation with JuMP.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.xlogx-Tuple{Matrix{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.xlogx","text":"xlogx(x::Matrix)\nreturn  x .* log.(x + eps() * (x .<= 0))\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.z_diag-Tuple{Vector{T} where T}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.z_diag","text":"z_diag(x::Vector)\n\nReturns a zero matrix with the input vector as diagonal.\n\n\n\n\n\n","category":"method"},{"location":"#Dupuy_Galichon_Matching.z_diag_f-Tuple{Matrix{T} where T, Function}","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Dupuy_Galichon_Matching.z_diag_f","text":"z_diag_f(x::Matrix, f::Function)\n\nMimick the matlab function diag(). \nApply a function on a Matrix columns (such as mean() or var()). \nReturns a zero matrix with the computed values in its diagonal.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","title":"Our Replication of Dupuy and Galichon \"Personality Traits and the Marriage Market\" (2014)","text":"end","category":"page"}]
}
