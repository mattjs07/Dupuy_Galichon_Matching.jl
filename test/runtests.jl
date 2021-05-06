using Test, Dupuy_Galichon_Matching, DataFrames, Statistics, LinearAlgebra

@test Dupuy_Galichon_Matching.convertion("0,01") == 0.01 

df = DataFrame(A = 1:4, B = ["1,01", "1.1", "2", "10,2"])
Dupuy_Galichon_Matching.convertion_df!(df)
@test  df == DataFrame(A = [1.0,2.0,3.0,4.0], B = [1.01, 1.1, 2, 10.2])


@test Dupuy_Galichon_Matching.z_diag_f([1:2 4:5 6:7], sum) == [3.0 0.0 0.0; 0.0 9.0 0.0; 0.0 0.0 13.0]