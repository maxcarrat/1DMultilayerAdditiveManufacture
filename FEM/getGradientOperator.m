function [ Bgl ] = getGradientOperator( problem, coords )
%GETGRADIENTOPERATOR get the gradient operator of the global problem

    %element gradient operator
    Bel = @(x) problem.basis_fun(x, 1, 1);
    
    %global gradient operator 
    Bgl = zeros(problem.N, 1);
    
    for e=1:problem.N
        ldof = 2;
        X1 = problem.coords(e);
        X2 = problem.coords(e+1);
        
        Bgl(e) = Bgl(e) + problem.B_map(X1,X2) * Bel(coords);
    end

end

