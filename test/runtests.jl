using Test, Dupuy_Galichon_Matching

@test hello("Julia") == "Hello, Julia"
@test domath(2.0, 5.0) ≈ 7.0

@test convertion("0,01") == 0.01