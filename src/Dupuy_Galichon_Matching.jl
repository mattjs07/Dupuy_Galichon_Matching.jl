module Dupuy_Galichon_Matching 

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using LinearAlgebra
using JuMP
using Ipopt

"""
    mm(x::Matrix, y::Matrix)

Performs matrice multiplication using loops. \\
Usefull for optimisation with JuMP.
"""
function mm(x::Matrix, y::Matrix)

    rowsx, colsx = size(x)
    rowsy, colsy = size(y)
    theMatrixProduct = zeros(rowsx, colsy)

for row in 1:rowsx
  for col in 1:colsy
    theSum = 0
    for k in 1:colsx
      theSum = theSum + x[row, k] * y[k, col]
    end
    theMatrixProduct[row, col] = theSum
    end
end
return theMatrixProduct
end


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
        x[!,i] = convertion.(x[!,i])
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


"""
    CovVarEdsi(X::Matrix,Y::Matrix )

Returns Cov _XY _XY , Cov _XX _XX, Cov _YY _YY, Cov _XX _YY
"""
function CovVarEdsi(X::Matrix,Y::Matrix )
    if size(X) != size(Y)
        error("dimensions of X != dimensions of Y")
    else
        N, Kx = size(X)
        N, Ky = size(Y)

        phi_xy = Array{Any,2}(undef,N,0)
        for i = 1:Kx
            for j = 1:Ky
             
              k = (i-1)*Ky + j
              phi_xy = hcat(phi_xy, X[:,i].*Y[:,j])
                      
            end
        end

        phi_xx = Array{Any,2}(undef,N,0)
        for i = 1:Kx
                for j = 1:Kx
                      
                      k = (i-1)*Kx + j
                      if i != j 
                        xx = 0*(X[:,i].*X[:,j])
                      else
                        xx = X[:,i].*X[:,j]
                      
                      end
                      phi_xx = hcat(phi_xx, xx)             
                end
        end

        phi_yy = Array{Any,2}(undef,N,0)
        for i = 1:Ky
                for j = 1:Ky
                          
                      k = (i-1)*Ky + j
                      if i != j 
                        yy = 0*(Y[:,i].*Y[:,j])
                      else
                        yy = Y[:,i].*Y[:,j]
                      
                      end
                      phi_yy = hcat(phi_yy, yy)             
                end
        end
        
        Cov_XY_XY = zeros(Kx*Ky,Kx*Ky)
        for k = 1:(Kx*Ky)
            for l = 1:(Kx*Ky)
                  
                Cov_XY_XY[k,l] = phi_xy[:,k]'*phi_xy[:,l]/N
                        
             end
        end

        Cov_XX_XX = zeros(Kx*Kx,Kx*Kx)
        for k = 1:(Kx*Kx)
            for l = 1:(Kx*Kx)
                  
                Cov_XX_XX[k,l] = phi_xx[:,k]'*phi_xx[:,l]/N
                        
             end
        end

        Cov_YY_YY = zeros(Ky*Ky,Ky*Ky)
        for k = 1:(Ky*Ky)
            for l = 1:(Ky*Ky)
                  
                Cov_YY_YY[k,l] = phi_yy[:,k]'*phi_yy[:,l]/N
                        
             end
        end

        Cov_XX_YY = zeros(Kx*Kx,Ky*Ky)
        for k = 1:(Kx*Kx)
            for l = 1:(Ky*Ky)
                  
                Cov_XX_YY[k,l] = phi_xx[:,k]'*phi_yy[:,l]/N
                        
             end
        end
    
    return Dict(:Cov_XY_XY => Cov_XY_XY , :Cov_XX_XX => Cov_XX_XX, :Cov_YY_YY => Cov_YY_YY, :Cov_XX_YY => Cov_XX_YY)
    
    end
        
        
end

"""
    EmpCov(X::Matrix,Y::Matrix)

Returns the empirical covariance of two matrices.
"""
function EmpCov(X::Matrix,Y::Matrix)

    N,M = size(X)
    
    Cov = zeros(M, M)
    
    for i = 1:N
         Cov = Cov + X[i,:] .* Y[i,:]'
    end
    
    return Cov/N
end

"""
    length_(X::Matrix)

Mimick matlab length() function by returning the length of the biggest dimension.
"""
function length_(X::Matrix)
    a,b = size(X)
    return max(a,b)

end

"""
    discretePdf(X::Matrix, Y::Matrix)

Returns two vectors of length equal to the matrices highest dimension (N). \\
Each element equals to 1 divided by the number of rows.
Used in OptimalPi() function.
"""
function discretePdf(X::Matrix, Y::Matrix)
    N = length_(X)

    pdfX = 1/ length(X[:,1]) * ones(N,1)
    pdfY = 1/ length(Y[:,1]) * ones(N,1)

    return [pdfX, pdfY]
end

"""
    StartingValues(X::Matrix, Y::Matrix)  
    
    N,M = size(X)
    u0 = log(N) * ones(N,1)
    v0 = -ones(N,1)
    return u0, v0

Used in OptimalPi() function.
"""
function StartingValues(X::Matrix, Y::Matrix)

    N,M = size(X)
    u0 = log(N) * ones(N,1)
    v0 = -ones(N,1)
    return u0, v0
end

"""
    OptimalPi(auxVar::Matrix, X::Matrix, Y::Matrix, T::Int64)

Returns the Optimal Pi, used in ObjectiveFunction() for optimization.
"""
function OptimalPi(auxVar::Matrix, X::Matrix, Y::Matrix, T::Int)

    maxIter = 10
    pdfX,pdfY = discretePdf(X,Y);
    
    r = exp.(auxVar/ T)
    
    normConstR = sum(sum(r,dims = 1))  
    r = r/normConstR
    
    u0, v0 = StartingValues(X,Y)
    a0 = exp.(u0)
    b0 = exp.(v0)
    
    ak = a0
    bk = b0
    
    for j in 1:maxIter
        bk = pdfY./(r'*ak)
        ak = pdfX./(r*bk)
    end
    
    normConst = ak' * r * bk
    
    Pi  = r .* (ak * bk')
    Pi  = Pi ./ normConst
    
    return Pi
    end
    
    """
        xlogx(x::Matrix)
        return  x .* log.(x + eps() * (x .<= 0))
    """
    function xlogx(x::Matrix)
    
        res = x .* log.(x + eps() * (x .<= 0))   ##
    end
    
    """
        ObjectiveFunction(X::Matrix,Y::Matrix,sigmaHat::Matrix,T::Int64,A... )
    Objective function used to compute the optimal Affinity Matrix. \\
    A is splat as it is the argument over which we arg min using JuMP.
    """
    function ObjectiveFunction(X::Matrix,Y::Matrix,sigmaHat::Matrix,T::Int64, colX::Int, A... )
    
        A = reshape(collect(A), (colX, colX)) #colX is the number of columns of X (same as Y)
        auxVar = X * A * Y'
    
        Pi = OptimalPi(auxVar,X,Y,T)
    
        TwistedTrace = sum(sum(auxVar .* Pi))
        entropy = sum(sum(-xlogx(Pi), dims = 1))
        OptCov = X' * Pi * Y
    
        f = TwistedTrace + T * entropy - tr(sigmaHat' *  A)
        G = OptCov - sigmaHat;
    
        return f   # Might want to get rid go G as the userdefined function to be optimized must return a Scalar !
    end


"""
AffinityMatrix(df::DataFrame)
Function allowing to compute the Affinity Matrix between two groups X and Y. \\
For m characteristics per group, returns a m x m Affinity Matrix. \\
!!! The input Dataframe contains columns of characteristics which must be ordered as follow : \\
> The Number of characteristics for X must be equal to the number of characteristics of Y. \\
> If the total number of columns = n, then columns 1:n/2 are X's and columns n/2+1:n are Y's. \\
> In each group the column must be ordered alike. \\
> e.g A Dataframe with 3 characteristics, the columns are arranged as :  [X.education, X.age, X.sex , Y.education, Y.age, Y.sex ] \\
"""
function AffinityMatrix(df::DataFrame)
convertion_df!(df)
col = size(df)[2]
m =  convert(Int, col/2)
X = select(df, 1:m )
Y = select(df, (m+1):col)
namesX = names(X)
namesY = names(Y)
X= Matrix(X)
Y= Matrix(Y)
zX = z_diag_f(X, mean)
X = X - ones(size(X)[1],size(X)[2]) * zX 
zY = z_diag_f(Y, mean)
Y = Y - ones(size(Y)[1],size(Y)[2]) * zY 
covX = cov(X)
sdX = sqrt.(covX[diagind(covX)])
SscaleX = z_diag(sdX)
covY = cov(Y)
sdY = sqrt.(covY[diagind(covY)])
SscaleY = z_diag(sdY)
X = X*inv(SscaleX); # divide each entry by the SD of its variable ==> standardizing
Y = Y*inv(SscaleY);
roX , colX = size(X);
roY , colY = size(Y);
CVE = CovVarEdsi(X,Y)
Cov_XY_XY, Cov_XX_XX, Cov_YY_YY, Cov_XX_YY =  CVE[:Cov_XY_XY], CVE[:Cov_XX_XX], CVE[:Cov_YY_YY], CVE[:Cov_XX_YY]
sigmaHat = EmpCov(X,Y)
roSig , colSig = size(sigmaHat);
T = 1
maxIter = 2000
m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter" => maxIter, "tol" => 1e-08    ))
set_optimizer_attribute(m, MOI.Silent(), true)
@variable(m, A[1:colX,1:colX]) 
JuMP.register(m, :ObjectiveFunction, nv, ObjectiveFunction, autodiff=true)
@NLobjective(m, Min, ObjectiveFunction(X, Y, sigmaHat, T, A...) ) #Does not work here 
@show JuMP.optimize!(m)
Aopt = @show value(A)
Aopt = NamedArray(A, (namesX, namesY))
return Aopt
end

end