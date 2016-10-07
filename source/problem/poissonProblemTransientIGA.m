%% Poisson Problem 1D transient IGA
function problem = poissonProblemTransientIGA(coords, rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, neumannBoundaryConditionValue, k, heatCapacity, time, ...
    knotVector, p, refinementDepth)

            % number of Elements
            N = size(knotVector, 2) - 2*p - 1;

            % LM = location matrix
            % it maps the shape functions local to each element to a global unknown
            % index
            LM = locationMapIGA(N, p);

            B_map = @(X1, X2) 2/(X2 - X1);
            F_map = @(X1, X2) (X2 - X1)/2;
            
            % Dirichlet BCs
            dirichlet_bc = [];
            
            if size(leftDirichletBoundaryConditionValue(time))~=0
                dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
            end
            if size(rightDirichletBoundaryConditionValue(time))~=0
                dirichlet_bc = [dirichlet_bc; LM(N,end) rightDirichletBoundaryConditionValue(time)];
            end
            
            % Neumann BCs
            if numel(neumannBoundaryConditionValue) ~= 0
                neumann_bc = [LM(N,2) neumannBoundaryConditionValue];
            end
            
            gdof = max(max(LM));
            
            penalty = 1.0e+12;

            problem = struct('LM', LM,  'B_map', B_map, 'F_map', F_map, 'refinementDepth', refinementDepth, ...
                'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', rhs, 'coords', coords, 'penalty', penalty, ...
                'k', k, 'heatCapacity', heatCapacity, 'time', time, 'neumann_bc', neumann_bc, 'knotVector', knotVector,...
                'p', p );
end 


