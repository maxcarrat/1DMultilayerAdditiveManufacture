function [poissonTransientProblemEnriched] = poissonNonLinearProblemTransientPODEnrichment(coords, rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, k, heatCapacity, time, refinementDepth, reductionOperator)
%POISSONNONLINEARPROBLEMTRANSIENTPODENRICHMENT struct of a Thermal transient non-linear
%problem using enriched POD basis functions

%% Standard problem

% number of Elements
N = size(coords, 2) - 1;

% basis function and location map
[basis_fun, LM] = locationMap(N);

%non-linear bilinear opertaor
B = @(u, v, T, X1, X2)  integral( @(xi) k(T(mapLocalToGlobal(xi, X1, X2),0)) .* u(xi, 1.0)...
                .* v(xi, 1.0), -1, 1);
            
%rhs
F = @(v, X1, X2) integral( @(xi) rhs(mapLocalToGlobal(xi, X1, X2), time)...
                .* v(xi, 0.0), -1, 1);
            
%capacity term operator
M = @(u, v)  integral( @(xi) heatCapacity .* u(xi, 0.0)...
                .* v(xi, 0.0), -1, 1);
            
%forward and backward 1D map
B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

%Dirichlet BCs 
dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
dirichlet_bc = [dirichlet_bc; LM(N,2) rightDirichletBoundaryConditionValue(time)];


%% Reduced Basis problem

%number of Modes
modes = size(reductionOperator, 2);

% basis function and location map
[rbBasis_fun, rbLM] = rbLocationMap(modes, 2);

%non-linear bilinear opertaor
rbB = @(u, v, T, X1, X2)  integral( @(xi) k(T(mapLocalToGlobal(xi, X1, X2),0)) .* u(xi)...
                .* v(xi), -1, 1);

%rhs
rbF = @(v, X1, X2) integral( @(xi) rhs(mapLocalToGlobal(xi, X1, X2), time)...
                .* v(xi), -1, 1);
            
%capacity term operator
rbM = @(u, v)  integral( @(xi) heatCapacity .* u(xi)...
                .* v(xi), -1, 1);
            
%Dirichlet BCs on the enriched DOFs
dirichlet_bc = [dirichlet_bc; 1 0.0];
dirichlet_bc = [dirichlet_bc; rbLM(2,2) 0.0];

gdof = max(max(rbLM));
cdof = max(max(LM));
edof = max(max(rbLM));

penalty = 1.0e+12;

%% Generate problem struct
poissonTransientProblemEnriched = struct('LM', LM, 'basis_fun', basis_fun,'rbBasis_fun', rbBasis_fun, 'B', B, 'B_map', B_map, 'F', F,...
    'F_map', F_map, 'M', M, 'dirichlet_bc', dirichlet_bc, 'refinementDepth', refinementDepth,...
    'N', N, 'gdof', gdof, 'cdof', cdof, 'edof', edof,...
    'rhs', rhs, 'penalty', penalty, 'rbB', rbB, 'rbM', rbM, 'rbF', rbF, 'coords', coords, 'k', k, 'heatCapacity',...
    heatCapacity, 'time', time, 'modes', modes, 'rbLM', rbLM, 'reductionOperator', reductionOperator);
end


