function [M, K, f] = assemblyMultiPhaseXIGASystem(problem, time, integrationOrder, integrationModalOrder,...
    solutionCoefficients, oldSolutionCoefficients)
%ASSEMBLYMULTIPHASEXIGASYSTEM: assembles the mass and the conductivity matrix and load vector
%   problem = definition of the boundary value problem
%   time = current time
%   integrationOrder = number of integration points

%% Allocate matrices

%global conductivity matrix
K = zeros(problem.gdof,problem.gdof);
%global capacity matrix
M = zeros(problem.gdof,problem.gdof);
%global load vector
f = zeros(problem.gdof, 1);

% IGA block
%IGA conductivity matrix
K_IGA = zeros(problem.IGAdof,problem.IGAdof);
%IGA capacity matrix
M_IGA = zeros(problem.IGAdof,problem.IGAdof);
%IGA load vector
f_IGA = zeros(problem.IGAdof, 1);

% XIGA block
%XIGA conductivity matrix
K_XIGA = zeros(problem.XIGAdof,problem.XIGAdof);
%XIGA capacity matrix
M_XIGA = zeros(problem.XIGAdof,problem.XIGAdof);
%XIGA load vector
f_XIGA = zeros(problem.XIGAdof, 1);

% Coupling block
%Coupling conductivity matrix
K_Coupling = zeros(problem.XIGAdof,problem.IGAdof);
%Coupling capacity matrix
M_Coupling = zeros(problem.XIGAdof,problem.IGAdof);


%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
[rGPXIGA, wGPXIGA] = gaussPoints( integrationModalOrder );

numberOfIntegrationPoints = length(rGP);
numberOfModalIntegrationPoints = length(rGPXIGA);

modes = problem.modes;
incrementModes = 0;

% On active elements use the refined domain as integration domain
refinedNodes = 2^problem.refinementDepth + 1;
integrationDomain = linspace(-1, 1, ceil(refinedNodes/problem.XN));
integrationCoefficients = linspace(-1, 1, ceil(refinedNodes/problem.XN));
subDomainShapeFunctionCoefficients = linspace(0, 1, ceil(refinedNodes/problem.XN));

for e=1:problem.N
    
    ldof = 2;
    
    Xi1 = problem.knotVector(e + problem.p);
    Xi2 = problem.knotVector(e + problem.p + 1);
    
    if e > (problem.N - problem.XN)
        
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            localCoords = mapParentToLocal(rGP(iGP), Xi1, Xi2);
            globalCoords = mapParentToGlobal(rGP(iGP), Xi1, Xi2, problem, e);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
            
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            inverseJacobianX_Xi = 1 / detJacobianX_Xi;
            
            %% Integrate IGA block
            %external heat source
            f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            %Capacity matrix
            M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGP(iGP), evaluateTemperature(e, localCoords, problem, solutionCoefficients), ...
                evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients))...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            %Diffusion matrix
            K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, evaluateTemperature(e, localCoords, problem, solutionCoefficients))...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGP(iGP) * inverseJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            elementEnrichedIndex = e - (problem.N - problem.XN);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = 2;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            modalDofs = length(indexLocalEnrichedNodes)*modes;
            
            for integrationSubDomainIndex =...
                    1 : ceil(refinedNodes/problem.XN)-1
                
                if rGPXIGA(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXIGA(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*floor(refinedNodes/problem.XN)+1:(elementEnrichedIndex-1)*floor(refinedNodes/problem.XN)...
                        + ceil(refinedNodes/problem.XN),:);
                    
                    [F, G] = PODModesAndDerivativesIGA( problem, rGPXIGA(iGP), modes, PODCoefficients,...
                        integrationCoefficients, integrationSubDomainIndex, indexLocalEnrichedNodes, e );
                    
                    globalGP = mapLocalToGlobal(rGPXIGA(iGP), Xi1, Xi2);
                    
                    %% Integrate Coupling block
                    %Capacity matrix
                    M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e))...
                        * F' * N * wGPXIGA(iGP) *  problem.F_map(Xi1,Xi2);
                    
                    %Diffusion matrix
                    K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) =...
                        K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e))...
                        * G' * B * wGPXIGA(iGP) * problem.B_map(Xi1,Xi2);
                    
                    %% Integrate XIGA block
                    %extrnal heat source
                    f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        problem.rhs(globalGP, time) * wGPXIGA(iGP) * problem.F_map(Xi1,Xi2);
                    
                    %Capacity matrix
                    M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e))...
                        * ( F' * F ) * wGPXIGA(iGP) * problem.F_map(Xi1,Xi2);
                    
                    %Diffusion matrix
                    K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationCoefficients, indexLocalEnrichedNodes, e))...
                        * ( G' * G ) * wGPXIGA(iGP) * problem.B_map(Xi1,Xi2);
                    break;
                end
            end
            
        end
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
    else
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            localCoords = mapParentToLocal(rGP(iGP), Xi1, Xi2);
            globalCoords = mapParentToGlobal(rGP(iGP), Xi1, Xi2, problem, e);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
            
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            inverseJacobianX_Xi = 1 / detJacobianX_Xi;
            
            %% Integrate IGA block
            %external heat source
            f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            %Capacity matrix
            M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGP(iGP), evaluateTemperature(e, localCoords, problem, solutionCoefficients), ...
                evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients))...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            %Diffusion matrix
            K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, evaluateTemperature(e, localCoords, problem, solutionCoefficients))...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGP(iGP) * inverseJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
        end
    end
    
end

K = [K_IGA, K_Coupling'; K_Coupling, K_XIGA];
M = [M_IGA, M_Coupling'; M_Coupling, M_XIGA];
f = [f_IGA; f_XIGA];


end


function [ projectedCoefficients ] = evaluateTemperature(e, x, problem, solutionCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(x), 2);
B = zeros(length(x), 2);

localCoords = x;

for k=1:length(x)
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
end

for i=1:length(x)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end

projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:));

end


function [ projectedCoefficients ] = evaluateTemperatureSubElements(e, x, problem, solutionCoefficients,...
    modes, subDomainShapeFunctionCoefficients, shapeFunctionCoefficients, indexLocalEnrichedNodes, element)
% EVALUATETEMPERATURESUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, length(solutionCoefficients(problem.N:end)) );

N = zeros(length(x), 2);

F = zeros(length(x), modes * numel(indexLocalEnrichedNodes));

localCoords = x;

for k=1:length(x)
    
    [N(k,:), ~] = shapeFunctionsAndDerivativesSubElements(localCoords(k), e,...
        subDomainShapeFunctionCoefficients);
    
    [F(k,:), ~] = PODModesAndDerivativesIGA(problem, localCoords(k), modes,...
        problem.reductionOperator, shapeFunctionCoefficients, e, indexLocalEnrichedNodes, element);
end

for i=1:length(localCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
end


projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);

end