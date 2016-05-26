function KE = rbLocalConductivityMatrix(problem, iMode, enrichedElementCoords)
%   KE = RBLOCALCONDUCTIVITYMATRIX(problem, refinementDepth) computes the conductivity local to one element
%   problem = definition of the boundary value problem
%   enrichedElementCoords = coordinates of the refined element

    ldof = 2;
    
    KE = zeros(ldof, ldof);
    
    for i=1:ldof
        for j=1:ldof
            KE(i,j) = problem.rbB(@(x)problem.rbBasis_fun(x, i, iMode, 1.0, problem, enrichedElementCoords ),...
                @(x)problem.rbBasis_fun(x, j, iMode, 1.0, problem, enrichedElementCoords));
        end
    end

end
