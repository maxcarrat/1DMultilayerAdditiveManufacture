%% Poisson Problem 1D transient XIGA
function problemXIGA = poissonProblemTransientXIGA(coords, activeNumberOfCPs,...
    rhs, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, neumannBoundaryConditionValue, k, heatCapacity, time, ...
    knotVector, p, refinementDepth, PODRefinementDepth, reductionOperator)

% number of Elements
numberOfElements = size(knotVector, 2) - 2*p - 1;
numberOfEnrichedElements = activeNumberOfCPs;

%number of Modes
modes = size(reductionOperator, 2);

% LM = location matrix
% it maps the shape functions local to each element to a global unknown
% index
LM = locationMapIGA(numberOfElements, p);
LMCoupling = locationMapIGACoupling(numberOfElements, numberOfEnrichedElements, p);
LMEnriched = locationMapXIGA(modes, numberOfEnrichedElements, p)';

B_map = @(X1, X2) 2/(X2-X1);
F_map = @(X1, X2) (X2-X1)/2;

dirichlet_bc = [];

%Dirichlet BCs of the strandard mesh
if size(leftDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; 1 leftDirichletBoundaryConditionValue(time)];
end
if size(rightDirichletBoundaryConditionValue(time))~=0
    dirichlet_bc = [dirichlet_bc; LM(numberOfElements,end) rightDirichletBoundaryConditionValue(time)];
end

% Neumann BCs
if numel(neumannBoundaryConditionValue) ~= 0
    neumann_bc = [LM(numberOfElements,end) neumannBoundaryConditionValue];
end

gdof = max(max(LM));
IGAdof = numberOfElements + 1;
XIGAdof = (2.^PODRefinementDepth) * modes;

penalty = 1.0e+15;

problemXIGA = struct('LM', LM, 'LMC', LMCoupling, 'LME', LMEnriched, 'B_map', B_map,...
    'F_map', F_map, 'dirichlet_bc', dirichlet_bc, 'neumann_bc', neumann_bc, 'refinementDepth', refinementDepth,...
    'N', numberOfElements, 'XN', numberOfEnrichedElements, 'gdof', gdof, 'IGAdof', IGAdof, 'XIGAdof', XIGAdof,...
    'rhs', rhs, 'penalty', penalty, 'coords', coords, 'k', k, 'heatCapacity', heatCapacity, 'knotVector', knotVector, ...
    'time', time, 'modes', modes,  'reductionOperator', reductionOperator, 'p', p );
end

