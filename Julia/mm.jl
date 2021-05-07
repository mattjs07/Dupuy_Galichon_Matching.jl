x = rand(5, 5)
y = rand(5, 5)


function mm(x::Matrix, y::Matrix)

    rowsx, colsx = size(x)
    rowsy, colsy = size(y)
    theMatrixProduct = zeros(rowsx, colsy)

for row in 1:rowsx
  for col in 1:colsy
    theSum = 0
    for k in 1:colsx
      theSum = theSum + x[row, k] * y[k, col]
    end
    theMatrixProduct[row, col] = theSum
    end
end
return theMatrixProduct
end

mm(x,y) .≈ x * y