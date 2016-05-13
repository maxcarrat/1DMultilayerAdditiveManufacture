function x = mapLocalToGlobal(xi, X1, X2)
    x = xi .* (X2-X1)/2 + (X2+X1)/2;
end