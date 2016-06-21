function M = rbLocalCapacityMatrix(problem, enrichedElementCoords)
%   M = RBLOCALCAPACITYMATRIX(problem) computes the load vector local to one element
%   problem = definition of the boundary value problem
%   enrichedElementCoords = coordinates of the refined element

    ldof = 2;
    numberOfModes = problem.modes;
    
    %sub-matrices
    M_FEM = zeros(ldof, ldof);
    M_enr = zeros(numberOfModes*ldof, numberOfModes*ldof);
    M_coupling = zeros(numberOfModes*ldof, ldof);
    
    for i=1:ldof
        for j=1:ldof
            M_FEM(i,j) = problem.rbM(@(x)problem.localBasis_fun(x, i, 1.0, problem, enrichedElementCoords),...
                @(x)problem.localBasis_fun(x, j, 1.0, problem, enrichedElementCoords));
        end
    end
    
    for i=1:ldof
        for iMode = 1:numberOfModes
            for j=1:ldof
                for jMode = 1:numberOfModes
                    M_enr((i-1)*numberOfModes + iMode,(j-1)*numberOfModes + jMode) = problem.rbM(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 0.0, problem, enrichedElementCoords ),...
                        @(x)problem.xFEMBasis_fun(x, j, jMode, 0.0, 0.0, problem, enrichedElementCoords));
                end
            end
        end
    end
    
    for j=1:ldof
        for i=1:ldof
            for iMode = 1:numberOfModes
                M_coupling((i-1)*numberOfModes + iMode,j) = problem.rbM(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 0.0, problem, enrichedElementCoords ),...
                    @(x)problem.localBasis_fun(x, j, 1.0, problem, enrichedElementCoords));
            end
        end
    end
    
    M = [M_FEM, M_coupling'; M_coupling, M_enr];

end
