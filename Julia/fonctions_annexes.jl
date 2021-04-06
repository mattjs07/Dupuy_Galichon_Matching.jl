using DataFrames
using DataFramesMeta


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

function convertion_df!(x::DataFrame)
    for i in 1:size(x)[2]
        x[i] = convertion.(x[i])
    end      
end

function z_diag_f(x::Matrix, f::Function)
    z = zeros(size(x)[2],size(x)[2]) 
    f_x = f.(eachcol(x))
    z[diagind(z)] = f_x
    return(z)
end

function z_diag(x::Vector)
    z = zeros(length(x),length(x))
    z[diagind(z)] = x
    return(z)
end

function CovVarEdsi(X::Matrix,Y::Matrix )
    if size(X) != size(Y)
        error("dimensions of X != dimensions of Y")
    else
        N = size(X)[1]
        K = size(X)[2]

        for i = 1:K
            for j = 1:K
             
              k = (i-1)*K + j;
              phi_xy(:,k) = X(:,i).*Y(:,j);
                         
            end
        end
        
    return phi_xy(:,k)

    















    end

end





