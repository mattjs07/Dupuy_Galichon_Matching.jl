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