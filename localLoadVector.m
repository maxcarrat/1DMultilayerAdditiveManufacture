function FE = localLoadVector(e, problem)
% FE = LOCALLOADVECTOR(coords, e, p, problem) computes the load vector local to one element
%   e = index of the element associated to the local load vector
%   problem = definition of the boundary value problem. 
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    ldof = 2;
    FE = zeros(ldof, 1);
    
    for i=1:ldof
        FE(i,1) = problem.F(@(x, d)problem.basis_fun(x, i, d), X1, X2);
    end


end