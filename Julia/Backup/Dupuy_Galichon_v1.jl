using CSV
using DataFrames
using DataFramesMeta
using Statistics
using LinearAlgebra

cd("C:/Users/matti/Desktop/M2_S2/Comp_Econ/Term_Project/Folder_1/Julia")
include("fonctions_annexes.jl")

df = DataFrame(CSV.File("BDG_w2040_27jun13.csv"))

deletecols!(df, [:gelukv, :gelukm,:Agem, :Agev])

convertion_df!(df)
df

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

CVE = CovVarEdsi(X,Y)
Cov_XY_XY, Cov_XX_XX, Cov_YY_YY, Cov_XX_YY =  CVE[:Cov_XY_XY], CVE[:Cov_XX_XX], CVE[:Cov_YY_YY], CVE[:Cov_XX_YY]

sigmaHat = EmpCov(X,Y)



# optimizes the objective function W(A) - trace(sigmaHat' * A) and returns the optimal matrix A
T = 1
maxIter = 2000
maxFunEvals = 10000

startMatrix = ones(colX, colX)


normConstraint = 100000;

# UNConstrained optimization, will have to use Optim.jl
#What is the objective function ? 
#There are parameters fed to the optimizing function, translate to Julia. 
# PRB :: they return the hessian (which Optim cannot), however JuMP can? Can we use JuMP for unconstrained optim? --> YES !! 
#However the documentation says it is not very good at it. And dunno if our ObjFund is automatically diferentiable ... 
# And in JuMP I know for sure we can feed Array variables 
# Still :::: What is the Godamn ObjectiveFunction ??? 
# Alternatives = Optim.jl (but cannot return hessian)
# Alternatives = NLopt.jl (cannot return hessian either ... )
# Hessian is a 100 x 100 matrix
using JuMP
using Ipopt

m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter" => maxIter, "tol" => 1e-08    ))
set_optimizer_attribute(m, MOI.Silent(), true)
@variable(m, A[1:colX,1:colX]) 
nv = colX^2
JuMP.register(m, :ObjectiveFunction, nv, ObjectiveFunction, autodiff=true)
@NLobjective(m, Min, ObjectiveFunction(X,Y,sigmaHat,T,A...) ) #Does not work here 
@show JuMP.optimize!(m)


using Optim

opti = optimize(A -> ObjectiveFunction(X,Y,sigmaHat,T,A), startMatrix , NelderMead()) #cannot compute the gradient and Hessian because it is in matrix form ...
                                                                                            # In matrix differentation, you need "Tensors" ! 









# The first reported table "Affinity.tex" reports the result of the optimization (no need hessian here)
# However, need the Hessian for what comes next : used as the FIsher Info
# Can we compute Hessian outside optimization ? 
