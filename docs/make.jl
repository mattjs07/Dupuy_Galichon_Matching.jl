push!(LOAD_PATH,"../src/")
using Documenter, Dupuy_Galichon_Matching

makedocs(modules = [Dupuy_Galichon_Matching], sitename = "Dupuy_Galichon_Matching.jl")

deploydocs(repo = "github.com/mattjs07/Dupuy_Galichon_Matching.jl.git", devbranch = "main")
