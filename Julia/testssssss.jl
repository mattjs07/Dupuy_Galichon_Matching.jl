phi_xy = Array{Any,2}(undef,N,0)
        for i = 1:Kx
            for j = 1:Ky
             
              k = (i-1)*Ky + j
              phi_xy = hcat(phi_xy2, X[:,i].*Y[:,j])
                         
            end
        end

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




for i = 1:Kx     # What this does = creates a vector 1:100
    for j = 1:Ky

      k = (i-1)*Ky + j
      println(k)
    end
end


for i = 1:Kx     # What this does = creates a vector 1:100
    for j = 1:Ky

      k = (i-1)*Ky + j
      println(k)
    end
end


a, b, c = 1 , 2 ,3

a,b,c,d = CVE[:Cov_XY_XY], CVE[:Cov_XX_XX], CVE[:Cov_YY_YY], CVE[:Cov_XX_YY]





f(x...) = reshape(collect(x), (10, 10))

A = rand(10,10)
f(Z...)

reshape(Z, (5,20))

ObjectiveFunction(X,Y,sigmaHat,T,A...)


function ObjectiveFunction(X::Matrix,Y::Matrix,sigmaHat::Matrix,T::Int64,A::Matrix )  ## Not sure regarding the type of A

  A = reshape(collect(A), (colX, colX))
  # Need to redefine the input A, using splatting ! (see .jpg in /Julia)
  auxVar = X * A * Y'

  Pi = OptimalPi(auxVar,X,Y,T)

  TwistedTrace = sum(sum(auxVar .* Pi))
  entropy = sum(sum(-xlogx(Pi), dims = 1))
  OptCov = X' * Pi * Y

  f = TwistedTrace + T * entropy - tr(sigmaHat' * A)
  G = OptCov - sigmaHat;

  return f, G
end


model = Model(Ipopt.Optimizer)
@variable(model, x[1:5] >= 0)
f(a::Vector,x...) = sum(x[i] *a[i] for i in 1:length(x))
register(model, :f, 5, f; autodiff = true)
@NLobjective(model, Min, f(x...))
@show JuMP.optimize!(model)


a= [1,2,3,4,5]
Z = rand(5,5)
model = Model(Ipopt.Optimizer)
@variable(model, x[1:5] >= 0)
function f(Z::Matrix,a::Vector,x...)
  x = reshape(collect(x),(1,5)) 
  x * Z * a
end
register(model, :f, 5, f; autodiff = true)
@NLobjective(model, Min, f(x...))
@show JuMP.optimize!(model)

function Obj(A... )

  A = reshape(collect(A), (colX, colX))

end