function ME = rbLocalCapacityMatrix(problem, iMode, enrichedElementCoords)
%   ME = RBLOCALCAPACITYMATRIX(problem) computes the load vector local to one element
%   problem = definition of the boundary value problem
%   enrichedElementCoords = coordinates of the refined element

    ldof = 2;
    
    ME = zeros(ldof, ldof);
    
    for i=1:ldof
        for j=1:ldof
            ME(i,j) = problem.rbM(@(x)problem.rbBasis_fun(x, i, iMode, 0.0, problem, enrichedElementCoords),...
                @(x)problem.rbBasis_fun(x, j, iMode, 0.0, problem, enrichedElementCoords));
        end
    end

end
