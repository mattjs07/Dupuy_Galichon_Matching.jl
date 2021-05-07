using Test, Dupuy_Galichon_Matching, DataFrames, Statistics, LinearAlgebra

@test Dupuy_Galichon_Matching.convertion("0,01") == 0.01 

df = DataFrame(A = 1:4, B = ["1,01", "1.1", "2", "10,2"])
Dupuy_Galichon_Matching.convertion_df!(df)
@test  df == DataFrame(A = [1.0,2.0,3.0,4.0], B = [1.01, 1.1, 2, 10.2])


@test Dupuy_Galichon_Matching.z_diag_f([1:2 4:5 6:7], sum) == [3.0 0.0 0.0; 0.0 9.0 0.0; 0.0 0.0 13.0]

@test Dupuy_Galichon_Matching.z_diag(repeat([1],3)) == [ 1 0 0; 0 1 0; 0 0 1]

x = [ 1 2; 1 3]
y = [3 1; 2 1]
Cov_XX_XX = [1.0 0.0 0.0 6.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 6.5 0.0 0.0 48.5]
Cov_YY_YY =  [48.5  0.0  0.0  6.5; 0.0  0.0  0.0  0.0; 0.0  0.0  0.0  0.0; 6.5  0.0  0.0  1.0]
Cov_XY_XY =  [6.5  2.5  15.0   6.0; 2.5  1.0   6.0   2.5; 15.0  6.0  36.0  15.0; 6.0  2.5  15.0   6.5]
Cov_XX_YY =  [6.5  0.0  0.0  1.0; 0.0  0.0  0.0  0.0; 0.0  0.0  0.0  0.0; 36.0  0.0  0.0  6.5 ]
@test Dupuy_Galichon_Matching.CovVarEdsi(x,y) == Dict(:Cov_XY_XY => Cov_XY_XY , :Cov_XX_XX => Cov_XX_XX, :Cov_YY_YY => Cov_YY_YY, :Cov_XX_YY => Cov_XX_YY)

@test Dupuy_Galichon_Matching.EmpCov(x,y) ==  [2.5  1.0; 6.0  2.5]

z = [1 2 3 4 5; 1 2 3 4 5]
@test Dupuy_Galichon_Matching.length_(z) == 5

q = Dupuy_Galichon_Matching.discretePdf(x, z) 
a = Matrix{Float64}(undef,2,1)
a .= 0.5
@test q[1] == a
@test q[2] == a

s1, s2 =Dupuy_Galichon_Matching.StartingValues(x,z)
z1, z2 = Matrix{Float64}(undef,2,1), Matrix{Float64}(undef,2,1)
z1 .= 0.6931471805599453
z2 .= -1.0
@test s1 ≈ z1
@test s2 == z2

z = [ 1 2; 1 2]
x, y = z, z
t = 1
@test Dupuy_Galichon_Matching.OptimalPi(z,x,y,t) == [ 0.25  0.25; 0.25  0.25]

l = Matrix{Float64}(undef,2,2)
l .= 1.3862943611198906  
@test Dupuy_Galichon_Matching.xlogx([ 2 2; 2 2]) ≈ l

a = [ 1 2 3 4]
colX = 2

@test Dupuy_Galichon_Matching.ObjectiveFunction(x,y,z,t,colX,a...) ≈ 11.38629436111989


