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
            M_FEM(i,j) = problem.rbM(@(x)problem.basis_fun(x, i, 0.0),...
                @(x)problem.basis_fun(x, j, 0.0));
        end
    end
    
    for iMode = 1:numberOfModes
        for i=1:ldof
            for jMode = 1:numberOfModes
                for j=1:ldof
                    M_enr((iMode-1)*ldof + i,(jMode-1)*ldof + j) = problem.rbM(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 0.0, problem, enrichedElementCoords ),...
                        @(x)problem.xFEMBasis_fun(x, j, jMode, 0.0, 0.0, problem, enrichedElementCoords));
                end
            end
        end
    end
    
    for j=1:ldof
        for i=1:ldof
            for iMode = 1:numberOfModes
                M_coupling((iMode-1)*ldof + i,j) = problem.rbM(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 0.0, problem, enrichedElementCoords ),...
                    @(x)problem.basis_fun(x, j, 0.0));
            end
        end
    end
    
    M = [M_FEM, M_coupling'; M_coupling, M_enr];

end
