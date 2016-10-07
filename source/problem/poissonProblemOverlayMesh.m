function problem  = poissonProblemOverlayMesh( coords, activeNumberOfControlPoints, rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, neumannBoundaryConditionValue, p, k, heatCapacity, time )
%POISSONPROBLEMOVERLAYMESH Create a struct for a 1D Poisson problem on the
%overlay FEM mesh in multiscale h-IGA
%   Detailed explanation goes here

% number of Elements
N = size(coords, 2) - 1;

% LM = location matrix
% LCoupling = location matrix for coupling base and overlay mesh
% it maps the shape functions local to each element to a global unknown
% index
[basis_fun, LM] = locationMap(N);

LCoupling = locationMapOverlayCoupling(N, activeNumberOfControlPoints, p)

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

% Dirichlet BCs
dirichlet_bc = [];

if size(leftDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
end
if size(rightDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; LM(N,2) rightDirichletBoundaryConditionValue(time)];
end

% Neumann BCs
if numel(neumannBoundaryConditionValue) ~= 0
    neumann_bc = [LM(N,2) neumannBoundaryConditionValue];
end

gdof = max(max(LM));

penalty = 1.0e+12;

problem = struct('LM', LM, 'LCoupling', LCoupling, 'basis_fun', basis_fun,  'B_map', B_map, 'F_map', F_map, ...
    'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', rhs, 'coords', coords, 'penalty', penalty, ...
    'k', k, 'heatCapacity', heatCapacity, 'time', time, 'neumann_bc', neumann_bc);

end

