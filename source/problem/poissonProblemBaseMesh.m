%% Poisson Problem 1D transient IGA
function problem = poissonProblemBaseMesh(coords, rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, neumannBoundaryCondition, k, steelThermalConductivityDerivative,...
    heatCapacity, heatCapacityDerivative, time, ...
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
            neumann_bc = [LM(N,end) rhs(coords(end), time)];
            
            gdof = max(max(LM));
            
            penalty = 1.0e+15;

            problem = struct('LM', LM,  'B_map', B_map, 'F_map', F_map, 'refinementDepth', refinementDepth, ...
                'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', rhs, 'coords', coords, 'penalty', penalty, ...
                'k', k, 'kDerivative',steelThermalConductivityDerivative,  'heatCapacity', heatCapacity,...
                'heatCapacityDerivative', heatCapacityDerivative, 'time', time, 'neumann_bc',...
                neumann_bc, 'knotVector', knotVector, 'p', p );
end 


