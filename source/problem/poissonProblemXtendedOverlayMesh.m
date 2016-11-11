function problemXOverlay = poissonProblemXtendedOverlayMesh(coords,...
    numberOfEnrichmentDepth, rhs, neumann_bc, p, k, steelThermalConductivityDerivative,...
    heatCapacity, heatCapacityDerivative, reductionOperator, time)
%POISSONPROBLEMXTENDEDOVERLAYMESH Create a struct for a 1D Poisson problem on the
%overlay FEM mesh in multiscale POD X-IGA

% number of Elements
N = size(coords, 2) - 1;

%number of Modes
modes = size(reductionOperator, 2);

numberOfEnrichedDofs = 2^(numberOfEnrichmentDepth);

% it maps the shape functions local to each element to a global unknown
% index
[~, LM] = locationMap(N);
LMBaseCoupling = locationMapXtendedBaseOverlayCoupling(N,...
    modes, numberOfEnrichedDofs, p+1);
LMCoupling = locationMapXtendedOverlayCouplingMesh(N, numberOfEnrichedDofs);
LMEnriched = locationMapXtendedOverlayMesh(modes, numberOfEnrichedDofs);

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

%Dirichlet BCs of the overlay mesh
dirichlet_bc = [dirichlet_bc; 1 0.0];
dirichlet_bc = [dirichlet_bc; LM(N,end) 0.0];

Xdirichlet_bc = [];

%Dirichlet BCs of the Xtended mesh
for i=1:modes*numberOfEnrichedDofs
    Xdirichlet_bc = [Xdirichlet_bc; i 0.0];
end

% Neumann BCs
neumann_bc = [LM(N,end) rhs(coords(end), time)];


overlayDofs = N + 1;
Xdof = (2*numberOfEnrichedDofs-1) * modes;
gdof = overlayDofs + Xdof;

penaltyL2 = 1.0e+15;
penalty = 1.0e+12;

problemXOverlay = struct('LM', LM, 'LMC', LMCoupling, 'LME', LMEnriched, 'B_map', B_map,...
    'F_map', F_map, 'dirichlet_bc', dirichlet_bc, 'neumann_bc', neumann_bc,...
    'XN', numberOfEnrichedDofs, 'N', N, 'gdof', gdof, 'LMBC', LMBaseCoupling, 'Xdof', Xdof,...
    'overDofs', overlayDofs, 'rhs', rhs, 'penalty', penalty,'penaltyL2', penaltyL2, 'coords', coords, 'k', k, 'kDerivative', steelThermalConductivityDerivative,...
    'heatCapacity', heatCapacity, 'heatCapacityDerivative', heatCapacityDerivative, ...
    'time', time, 'modes', modes,  'reductionOperator', reductionOperator,'Xdirichlet_bc', Xdirichlet_bc );

end

