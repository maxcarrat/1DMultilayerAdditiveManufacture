function refinedLocalProblem = poissonProblemXFEM( coords, activeMeshSize, rhs, leftDirichletBoundaryConditionValue,...
                                rightDirichletBoundaryConditionValue, k, heatCapacity, time, refinementDepth, reductionOperator )
%POISSONPROBLEMTRANSIENTENRICHED The right boundary element of the global coarse problem
%is enriched using the POD modes from the training phase.

% number of Elements
numberOfElements = size(coords, 2) - 1;
numberOfEnrichedElements = activeMeshSize;

%number of Modes
modes = size(reductionOperator, 2);

[~, LM] = locationMap(numberOfElements);
LMCoupling = locationMapCoupling(numberOfElements, numberOfEnrichedElements);
LMEnriched = locationMapEnriched(modes, numberOfEnrichedElements);

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

%Dirichlet BCs of the strandard mesh
if size(leftDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
end
if size(rightDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; LM(numberOfElements,2) rightDirichletBoundaryConditionValue(time)];
end

%Dirichlet BCs on the enriched DOFs
dirichlet_bc = [dirichlet_bc; 1 0.0];
dirichlet_bc = [dirichlet_bc; LM(2,2) 0.0];


gdof = max(max(LM));
FEMdof = numberOfElements + 1;
XFEMdof = 2 * modes;

penalty = 1.0e+20;

refinedLocalProblem = struct('LM', LM, 'LMC', LMCoupling, 'LME', LMEnriched, 'B_map', B_map,...
    'F_map', F_map, 'dirichlet_bc', dirichlet_bc, 'refinementDepth', refinementDepth,...
    'N', numberOfElements, 'XN', numberOfEnrichedElements, 'gdof', gdof, 'FEMdof', FEMdof, 'XFEMdof', XFEMdof,...
    'rhs', rhs, 'penalty', penalty, 'coords', coords, 'k', k, 'heatCapacity',...
    heatCapacity, 'time', time, 'modes', modes,  'reductionOperator', reductionOperator);
end

