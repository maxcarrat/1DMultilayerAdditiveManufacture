function [ xFEMProblem ] = poissonProblemArmonicXFEM( coords, rhs, leftDirichletBoundaryConditionValue,...
                                rightDirichletBoundaryConditionValue, k, heatCapacity, time, numberOfModes )
%POISSONPROBLEMARMONICXFEM The right boundary element of the global coarse problem
%is enriched using the POD modes from the training phase.

% number of Elements
N = size(coords, 2) - 1;

%number of Modes

[basis_fun, LM] = locationMap(N);
[rbBasis_fun, rbLM] = armonicLocationMap(N, numberOfModes);

xFEMBasis_fun = @xFEMBasis;

rbB = @(u, v)  integral( @(xi) k .* u(xi)...
                .* v(xi), -1, 1);

rbF = @(v, X1, X2) integral( @(xi) rhs(mapLocalToGlobal(xi, X1, X2), time)...
                .* v(xi), -1, 1);

rbM = @(u, v)  integral( @(xi) heatCapacity .* u(xi)...
                .* v(xi), -1, 1);
            
B = @(u, v)  integral( @(xi) k .* u(xi, 1.0)...
                .* v(xi, 1.0), -1, 1);

F = @(v, X1, X2) integral( @(xi) rhs(mapLocalToGlobal(xi, X1, X2), time)...
                .* v(xi, 0.0), -1, 1);

M = @(u, v)  integral( @(xi) heatCapacity .* u(xi, 0.0)...
                .* v(xi, 0.0), -1, 1);

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

%Dirichlet BCs of the strandard mesh
if size(leftDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
end
if size(rightDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; LM(N,2) rightDirichletBoundaryConditionValue(time)];
end


gdof = N + 1 + N * (2*numberOfModes);
cdof = max(max(LM));
edof = max(max(rbLM));

penalty = 1.0e+10;

xFEMProblem = struct('LM', LM, 'basis_fun', basis_fun,'rbBasis_fun', rbBasis_fun, 'B', B, 'B_map', B_map, 'F', F,...
    'F_map', F_map, 'M', M, 'dirichlet_bc', dirichlet_bc,...
    'N', N, 'gdof', gdof, 'cdof', cdof, 'edof', edof, 'xFEMBasis_fun', xFEMBasis_fun,...
    'rhs', rhs, 'penalty', penalty, 'rbB', rbB, 'rbM', rbM, 'rbF', rbF, 'coords', coords, 'k', k, 'heatCapacity',...
    heatCapacity, 'time', time, 'modes', numberOfModes, 'rbLM', rbLM);
end