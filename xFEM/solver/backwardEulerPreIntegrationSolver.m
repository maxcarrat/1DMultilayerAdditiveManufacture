function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime]= backwardEulerPreIntegrationSolver(coords, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector,...
    refinementDepth, numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes, integrationOrder)
% BackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   rhs = external heat source
%   leftDirichletBoundaryConditionValue
%   rightDirichletBoundaryConditionValue
%   k = conductivity of the material
%   heatCapacity = heat capacity of the material
%   timeVector = vector of time steps for Backward Euler implicit scheme
%   refinementDepth = depth of the refined mesh
%   numberOfTrainingTimeSteps = number of time steps in the training phase

timeSteps=size(timeVector,2);

timeStepSize=max(timeVector)/( timeSteps );
temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

refinedTemperatureSolutions = zeros(2^refinementDepth+1, 1);
computationalTime = [];
localRefinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

previousMesh = coords;

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for layer = 1:numberOfTrainingLayers
    
    
    %Refinement
    refinedMesh = refineLayer(coords, refinementDepth, layer, numberOfLayers);
    %     refinedMesh = refineMesh(coords, refinementDepth, layer, numberOfTrainingLayers);
    
    %Project old solution onto the new mesh
    if layer > 1
        %         refinedTemperatureSolutions = L2projection(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature);
        refinedTemperatureSolutions = L2projectionLinearDistribution(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh,...
            previousMesh, initialTemperature, refinementDepth);
    else
        refinedTemperatureSolutions = zeros(size(refinedMesh, 2), 1);
        localRefinedTemperatureSolutions = zeros(size(refinedMesh, 2), timeSteps);
    end
    
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime);
        
        %Update and merge temperature into global domain
        refinedTemperatureSolutions = solveGlobalProblemGaussIntegration( refinedTemperatureSolutions, poissonTransientProblem,...
            currentTime, timeStepSize, integrationOrder );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
        
        localRefinedTemperatureSolutions(:, t+1) = getLayerSolution(refinedTemperatureSolutions, layer, numberOfLayers, coords);
        previousMesh = refinedMesh;
        
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
        %         internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
        
    end
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,3:numberOfTrainingLayers), numberOfPODModes);

%% Enriched mesh using RB

%Project old solution onto the new mesh
refinedDOFs = 2^refinementDepth;
temperatureSolutions = [refinedTemperatureSolutions(1:end-refinedDOFs);...
    refinedTemperatureSolutions(end); initialTemperature];

% %Generate Local problem
% [activeMesh, numberOfActiveElementsLayer] = getLayerActiveCoords(coords, layer, numberOfLayers);
%
% poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
%     rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);


for layer = (numberOfTrainingLayers+1):numberOfLayers
    
    [activeMesh, numberOfActiveElementsLayer] = getLayerActiveCoords(coords, layer, numberOfLayers);
    
    poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
    
    [K_Coupling, K_XFEM, M_Coupling, M_XFEM, F_XFEM, modalDofs] = integrateExtendedTerms(poissonTransientProblemEnriched, integrationOrder);
    
    %Project global solution onto the enriched modal space
    temperatureSolutions = eXtendedProjection(poissonTransientProblemEnriched,temperatureSolutions,...
        modes, initialTemperature);
    
    timeToGenerateAndSolveTheSystem = 0.0;
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * (t-1);
        
        %Generate Local problem
        poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
            rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
        
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        temperatureSolutions = solveEnrichedProblemPreIntegration( temperatureSolutions, poissonTransientProblemEnriched, currentTime,...
            timeStepSize, integrationOrder, K_Coupling, K_XFEM, M_Coupling, M_XFEM, F_XFEM, modalDofs );
        
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t) = postProcessingProjection(postProcessingCoords, poissonTransientProblemEnriched, temperatureSolutions,...
            modes, 0.0);
        heatFluxes(:, t) = postProcessingProjection(postProcessingCoords, poissonTransientProblemEnriched,temperatureSolutions,...
            modes, 1.0);
        %     [~, K, ~] = assembly(poissonTransientProblem);
        %     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

end

function [ projectedCoefficients ] = eXtendedProjection(problem, previousTemperature,...
    modes, initialTemperature)
% EXTENDEDPROJECTION project the previous solution onto the updated mesh at the
% new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   coarseMesh = initial/coarse mesh configuration
%   modes = number of enrichment modes
%   initialTemperature = initial temperature of the powder

[rGL, ~] = GaussLobatto(modes+1);
numberOfProjectionPoints = length(rGL);
projectedCoefficients = [];

for e=1:problem.N
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    if e > problem.N - problem.XN  % element is active
        
        if e == problem.N - problem.XN +1
            indexLocalEnrichedNodes = 2; %rhs node
        else
            indexLocalEnrichedNodes = [1, 2];
        end
        
        % rhs node
        projectedCoefficients = [projectedCoefficients; previousTemperature(e)];
        % lhs node
        projectedCoefficients = [projectedCoefficients; initialTemperature];
        
        % On active elements use the refined domain as integration domain
        integrationDomain = linspace(-1, 1, 2^problem.refinementDepth + 1);
        subDomainShapeFunctionCoefficients = linspace(0, 1, 2^problem.refinementDepth + 1);
        
        projectionOperator = [];
        projectionValues = [];
        
        for integrationSubDomainIndex=1:length(integrationDomain)-1
            
            Xi1 = integrationDomain(integrationSubDomainIndex);
            Xi2 = integrationDomain(integrationSubDomainIndex + 1);
            
            % Gauss-Lobatto point to project
            for iGL = 1:numberOfProjectionPoints
                
                if rGL(iGL) >= integrationDomain(integrationSubDomainIndex) && rGL(iGL) <= integrationDomain(integrationSubDomainIndex + 1)
                    
                    [N, ~] = shapeFunctionsAndDerivativesSubElements(rGL(iGL), integrationSubDomainIndex,...
                        subDomainShapeFunctionCoefficients);
                    [F, ~] = PODModesAndDerivatives(rGL(iGL), modes, problem.reductionOperator, subDomainShapeFunctionCoefficients,...
                        integrationSubDomainIndex,indexLocalEnrichedNodes);
                    
                    projectionOperator = [projectionOperator; N, F];
                    if iGL == 1
                        projectionValues = [projectionValues; 0.0];
                        
                    else
                        projectionValues = [projectionValues; 0.0];
                    end
                    
                end
            end
        end
        
        projectedCoefficientsEnriched = projectionOperator\projectionValues;
        projectedCoefficients = [projectedCoefficients; projectedCoefficientsEnriched(3:end)];
        
    else     % Element not active
        projectedCoefficients = [projectedCoefficients; previousTemperature(e)];
    end
end


end

function [ projectedCoefficients ] = postProcessingProjection(x, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTION project the previous solution onto the post-processing mesh.
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives


projectedCoefficients=zeros(size(x));
X1 = problem.coords(1);
X2 = problem.coords(2);

projectedCoefficients(x>=X1 & x<=X2) = postProcessingProjectionElement(1, x(x>=X1 & x<=X2), problem,... % first element
    solutionCoefficients, modes, derivative);

% loop over elements
for e=2:problem.N
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    projectedCoefficients(x>X1 & x<=X2) = postProcessingProjectionElement(e, x(x>X1 & x<=X2), problem,...
        solutionCoefficients, modes, derivative);
    
end

end

function [ projectedCoefficients ] = postProcessingProjectionElement(e, x, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTIONELEMENT project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);

X1 = problem.coords(e);
X2 = problem.coords(e+1);

if e > problem.N - problem.XN  % element is active
    
    if e == problem.N - problem.XN +1
        indexLocalEnrichedNodes = 2; %rhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    integrationDomain = linspace(-1, +1, 2^problem.refinementDepth + 1);
    subDomainShapeFunctionCoefficients = linspace(0, 1, 2^problem.refinementDepth + 1);
    
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    localCoordinates = mapGlobalToLocal( x, X1, X2);
    projectedCoefficients=zeros(size(localCoordinates));
    
    projectedCoefficients(localCoordinates>=Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(1,...
        localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), problem, integrationDomain,...
        solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients);
    
    for integrationSubDomainIndex=2:length(integrationDomain)-1
        
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
        
        projectedCoefficients(localCoordinates>Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(integrationSubDomainIndex,...
            localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), problem, integrationDomain,...
            solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients);
        
    end
    
else     % Element not enriched
    
    projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
    
    N = zeros(length(x), 2);
    B = zeros(length(x), 2);
    
    localCoords = mapGlobalToLocal( x, X1, X2);
    
    for k=1:length(x)
        [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    end
    
    if derivative == 0
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
    else
        for i=1:length(x)
            projectionOperator(i,1:size(B,2)) = B(i,:);
        end
    end
    
    projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:));
    
end


end


function [ projectedCoefficients ] = postProcessingProjectionSubElements(e, x, problem, integrationDomain, solutionCoefficients,...
    modes, derivative, subDomainShapeFunctionCoefficients)
% POSTPROCESSINGPROJECTIONSUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, length(solutionCoefficients(problem.N:end)) );

X1 = integrationDomain(e);
X2 = integrationDomain(e+1);
localCoords =  mapGlobalToLocal( x, X1, X2 );

N = zeros(length(x), 2);
B = zeros(length(x), 2);

F = zeros(length(x), modes);
G = zeros(length(x), modes);

for k=1:length(x)
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivativesSubElements(localCoords(k), e,...
        subDomainShapeFunctionCoefficients);
    [F(k,:), G(k,:)] = PODModesAndDerivatives(localCoords(k), modes, problem.reductionOperator,...
        subDomainShapeFunctionCoefficients, e, 2);
end

if derivative == 0
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
else
    for i=1:length(localCoords)
        projectionOperator(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
    end
end

projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);

end