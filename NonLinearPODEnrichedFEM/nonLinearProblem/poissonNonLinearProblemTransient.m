function [ nonLinearProblem ] = poissonNonLinearProblemTransient(coords, rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, k, heatCapacity, time)
%POISSONNONLINEARPROBLEMTRANSIENT struct of a Thermal transient non-linear
%problem

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
if size(leftDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
end
if size(rightDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; LM(N,2) rightDirichletBoundaryConditionValue(time)];
end

gdof = max(max(LM));

penalty = 1.0e+12;

nonLinearProblem = struct('LM', LM, 'basis_fun', basis_fun,'B', B, 'B_map', B_map, 'F', F,...
    'F_map', F_map, 'M', M, 'dirichlet_bc', dirichlet_bc, ...
    'N', N, 'gdof', gdof, 'rhs', rhs, 'penalty', penalty, 'coords', coords, 'k', k, 'heatCapacity',...
    heatCapacity, 'time', time);
end

