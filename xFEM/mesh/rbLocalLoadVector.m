function FE = rbLocalLoadVector(problem, iMode, enrichedElementCoords, e)
% FE = RBLOCALLOADVECTOR(problem, enrichedElementCoords, e) computes the load vector local to one element
%   problem = definition of the boundary value problem. 
%   enrichedElementCoords = coordinates of the refined element
%   e = index of the element associated to the local load vector

    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    ldof = problem.modes;
    FE = zeros(ldof, 1);
    
    for i=1:ldof
        FE(i,1) = problem.rbF(@(x)problem.rbBasis_fun(x, i, iMode, 0.0, problem, enrichedElementCoords), X1, X2);
    end


end