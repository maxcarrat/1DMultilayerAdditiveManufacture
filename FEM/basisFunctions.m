function y = basisFunctions(x, i, derivative)
    funs = {@getBasis, @getBasisDerivatives};
    
    y= funs{derivative+1}(x,i);
end



function y = getBasis(x, i)
   
    switch i
        case 1
            y = 0.5.*(1-x);

        case 2
            y = 0.5.*(1+x);
    end
end

function y = getBasisDerivatives(x, i)

    switch i
        case 1
            y = ones(size(x)).* -0.5;

        case 2
            y = ones(size(x)).* 0.5;
    end
end
