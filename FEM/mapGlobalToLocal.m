function Xi = mapGlobalToLocal(X, x1, x2)
    Xi = X .* 2/(x2-x1)  - (x2+x1) / (x2-x1);
end