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


[Cov_XY_XY, Cov_XX_XX, Cov_YY_YY, Cov_XX_YY] = CovVarEdsi(X,Y)

















