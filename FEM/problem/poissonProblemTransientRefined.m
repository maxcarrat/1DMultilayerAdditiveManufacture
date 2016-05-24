function refinedLocalProblem = poissonProblemTransientRefined( problemCoarse, coarseSolutions, refinementDepth )
%POISSONPROBLEMTRANSIENTREFINED The right boundary element of the global coarse problem
%is refined up to a user defined refinementDepth.The local problem is
%defined using as BCs the solutions on the coarse mesh.

localBoundaryCoords = [problemCoarse.coords(end-1) problemCoarse.coords(end)];

coordsRefined = linspace(localBoundaryCoords(1), localBoundaryCoords(2), 2^refinementDepth+1);

% number of Elements
N = size(coordsRefined, 2) - 1;

[basis_fun, LM] = locationMap(N);

B = @(u, v)  integral( @(xi) problemCoarse.k .* u(xi, 1.0) .* v(xi, 1.0), -1, 1);

F = @(v, X1, X2) integral( @(xi) problemCoarse.rhs(mapLocalToGlobal(xi, X1, X2), problemCoarse.time) .* v(xi, 0.0), -1, 1);

M = @(u, v)  integral( @(xi) problemCoarse.heatCapacity .* u(xi, 0.0) .* v(xi, 0.0), -1, 1);

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

dirichlet_bc = [dirichlet_bc; 1 coarseSolutions(end-1)];
dirichlet_bc = [dirichlet_bc; LM(N,2) problemCoarse.dirichlet_bc(2,2)];


gdof = max(max(LM));

penalty = 1.0e+12;

refinedLocalProblem = struct('LM', LM, 'basis_fun', basis_fun, 'B', B, 'B_map', B_map, 'F', F, 'F_map', F_map, 'M', M, ...
    'dirichlet_bc', dirichlet_bc, 'N', N, 'gdof', gdof, 'rhs', problemCoarse.rhs, 'coords', coordsRefined, 'penalty', penalty, ...
    'k', problemCoarse.k, 'heatCapacity', problemCoarse.heatCapacity, 'time', problemCoarse.time);
end

