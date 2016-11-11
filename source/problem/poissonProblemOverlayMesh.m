function problem  = poissonProblemOverlayMesh( coords, activeNumberOfControlPoints, rhs,...
    neumannBoundaryCondition, p, k, steelThermalConductivityDerivative,...
    heatCapacity, heatCapacityDerivative, time )
%POISSONPROBLEMOVERLAYMESH Create a struct for a 1D Poisson problem on the
%overlay FEM mesh in multiscale h-IGA

% number of Elements
N = size(coords, 2) - 1;

% LM = location matrix
% LCoupling = location matrix for coupling base and overlay mesh
% it maps the shape functions local to each element to a global unknown
% index
[basis_fun, LM] = locationMap(N);

% LCoupling = locationMapOverlayCoupling(N, activeNumberOfControlPoints, p);

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

% Dirichlet BCs
dirichlet_bc = [];

dirichlet_bc = [dirichlet_bc; 1 0.0];

dirichlet_bc = [dirichlet_bc; LM(N,2) 0.0];


% Neumann BCs
neumann_bc = [LM(N,end) rhs(coords(end), time)];


gdof = max(max(LM));
penaltyL2 = 1.0e+10;
penalty = 1.0e+08;

problem = struct('LM', LM, 'basis_fun', basis_fun,  'B_map', B_map, 'F_map', F_map, ...
    'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', rhs, 'coords', coords, 'penalty', penalty, ...
    'penaltyL2', penaltyL2, 'k', k, 'kDerivative', steelThermalConductivityDerivative, 'heatCapacity', heatCapacity,...
    'heatCapacityDerivative', heatCapacityDerivative, 'time', time, 'neumann_bc', neumann_bc);

end

