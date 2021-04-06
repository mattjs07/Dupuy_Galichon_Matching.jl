using CSV
using DataFrames
using DataFramesMeta
using Statistics

include("fonctions_annexes.jl")


cd("C:/Users/matti/Desktop/M2_S2/Comp_Econ/Term_Project/Folder_1/Julia")
df = DataFrame(CSV.File("BDG_w2040_27jun13.csv"))

#  df2 = select!(df, [Not(:gelukm, :gelukv, :Agem, :Agev)])
deletecols!(df, [:gelukv, :gelukm,:Agem, :Agev])

c = size(df)[2]
m =  convert(Int, c/2)
X = select(df, 1:m )
Y = select(df, (m+1):n)

mX = mean(X, dims=1)

mean.(eachcol(X))

convertion.(Float64, eachcol(X))

df[:,2] = convertion.(Float64,df[:,2])

df[:,2] = convert.(Float64, df[:,2])

for i in 1:n
    df[:,i] = convertion.(df[:,i])



end
df.consm = parse.(Float64, df[:,:])

convertion = function(x::Any) 
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




convertion(string("-10,33"))

parse(Float64, string(10))


typeof(string(10))

parse(Float64, 10 )

df

names(df)

convertion(df[1,1])
df[1,5] #le probleme c'est la virgule 

z = string("1,22")
occursin(",", z)

replace(z, ',', '.')

string(z[1],".",z[3:4])
z[2]
z[1] != ','
z[2] = '.'
z[1:1]

c = 1
zz = z[c]
while zz != ','
    c += 1
    zz = z[c]
end
c


string(x[1:(c-1)],".",x[(c+1):length(x)])



DecimalFormatSymbols symbols = new DecimalFormatSymbols();
symbols.setDecimalSeparator('.');
DecimalFormat format = new DecimalFormat("0.#");
format.setDecimalFormatSymbols(symbols);
float f = format.parse(str).floatValue();