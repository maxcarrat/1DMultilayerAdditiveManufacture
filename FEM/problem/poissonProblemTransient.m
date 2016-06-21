%% Poisson Problem 1D transient
function problem = poissonProblemTransient(coords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, time)

            % number of Elements
            N = size(coords, 2) - 1;

            % LM = location matrix
            % it maps the shape functions local to each element to a global unknown
            % index
            [basis_fun, LM] = locationMap(N);

            B = @(u, v)  integral( @(xi) k .* u(xi, 1.0) .* v(xi, 1.0), -1, 1);
            
            F = @(v, X1, X2) integral( @(xi) rhs(mapLocalToGlobal(xi, X1, X2), time) .* v(xi, 0.0), -1, 1);
                        
            M = @(u, v)  integral( @(xi) heatCapacity .* u(xi, 0.0) .* v(xi, 0.0), -1, 1);

            B_map = @(X1, X2) 2/(X2-X1);
            F_map = @(X1, X2) (X2-X1)/2;
            
            %Dirichlet BCs
            dirichlet_bc = [];
            
            if size(leftDirichletBoundaryConditionValue(time))~=0
                dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
            end
            if size(rightDirichletBoundaryConditionValue(time))~=0
                dirichlet_bc = [dirichlet_bc; LM(N,2) rightDirichletBoundaryConditionValue(time)];
            end
            
            gdof = max(max(LM));
            
            penalty = 1.0e+12;

            problem = struct('LM', LM, 'basis_fun', basis_fun, 'B', B, 'B_map', B_map, 'F', F, 'F_map', F_map, 'M', M, ...
                'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', rhs, 'coords', coords, 'penalty', penalty, ...
                'k', k, 'heatCapacity', heatCapacity, 'time', time);
end 


