function F = rbLocalLoadVector(problem, enrichedElementCoords, e)
% F = RBLOCALLOADVECTOR(problem, enrichedElementCoords, e) computes the load vector local to one element
%   problem = definition of the boundary value problem. 
%   enrichedElementCoords = coordinates of the refined element
%   e = index of the element associated to the local load vector

    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    ldof = 2;
    numberOfModes = problem.modes;
    
    %sub-matrices
    F_FEM = zeros(ldof, 1);
    F_enr = zeros(numberOfModes*ldof, 1);
    
    for i=1:ldof
            F_FEM(i,1) = problem.rbF(@(x)problem.basis_fun(x, i, 0.0), X1, X2);
    end
    
    for iMode = 1:numberOfModes
        for i=1:ldof
            F_enr((iMode-1)*ldof + i, 1) =...
                problem.rbF(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 0.0, problem, enrichedElementCoords), X1, X2);
        end
    end

    F = [F_FEM; F_enr];

end