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

refinedControlPoints = 2^problem.refinementDepth + problem.p;
% refinedControlPoints = 2^problem.refinementDepth + 1;
integrationDomain = linspace(-1, 1, ceil(refinedControlPoints/problem.XN));
subDomainShapeFunctionCoefficients = linspace(0, 1, ceil(refinedControlPoints/problem.XN));
subDomainKnotVector =  getOpenKnotVector( ceil(refinedControlPoints/problem.XN) , problem.p );

for e=1:problem.N
    
    ldof = problem.p + 1;
    
    Xp1 = problem.knotVector(e + problem.p);
    Xp2 = problem.knotVector(e + problem.p + 1);
    
    if e > (problem.N - problem.XN)
        
        
        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            localCoords = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
            globalCoords = mapParentToGlobal(rGPXIGA(iGP), Xp1, Xp2, problem, e);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
            
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            inverseJacobianX_Xi = 1 / detJacobianX_Xi;
            
            %% Integrate IGA block
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGPXIGA(iGP), evaluateTemperature(e, localCoords, problem, solutionCoefficients), ...
                evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients))...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            %Diffusion matrix
            K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, evaluateTemperature(e, localCoords, problem, solutionCoefficients))...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * inverseJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            elementEnrichedIndex = e - (problem.N - problem.XN);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = problem.IGAdof;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            modalDofs = length(indexLocalEnrichedNodes)*modes;
            
            for integrationSubDomainIndex =...
                    1 : ceil(refinedControlPoints/problem.XN)-1
                
                if rGPXIGA(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXIGA(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*floor(refinedControlPoints/problem.XN)+1:(elementEnrichedIndex-1)*floor(refinedControlPoints/problem.XN)...
                        + ceil(refinedControlPoints/problem.XN),:);
                    
                    x = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
                    
%                     Xsub1 = subDomainKnotVector(integrationSubDomainIndex);
%                     Xsub2 = subDomainKnotVector(integrationSubDomainIndex + 1);
%            
%                     xSub = mapParentToLocal(rGPXIGA(iGP), Xsub1, Xsub2);
% 
%                     [N, B] = shapeFunctionsAndDerivativesSubElements(rGPXIGA(iGP), integrationSubDomainIndex, N);
%           
% 
%                     [N, B] = BsplinesShapeFunctionsAndDerivatives(x, problem.p, subDomainKnotVector);

                    [N, B] = shapeFunctionsAndDerivativesIGASubElements( rGPXIGA(iGP), integrationSubDomainIndex,...
                        subDomainShapeFunctionCoefficients,  Xp1, Xp2,...
                        integrationDomain(integrationSubDomainIndex), integrationDomain(integrationSubDomainIndex+1), problem );
                    
                    [F, G] = PODModesAndDerivativesIGA( problem, x, modes, PODCoefficients,...
                        integrationDomain, integrationSubDomainIndex, indexLocalEnrichedNodes, e );
                  
%                     [F, G] = PODModesAndDerivativesGaussIntegration(rGPXIGA(iGP), modes, PODCoefficients,...
%                         integrationCoefficients, integrationSubDomainIndex, 2 );                    
                   
                    globalGP = mapParentToGlobal(rGPXIGA(iGP), Xp1, Xp2, problem, e);
                    
                    %% Integrate Coupling block
                    %Capacity matrix
                    M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain))...
                        * F' * N(end-problem.p:end) * wGPXIGA(iGP) *  detJacobianX_Xi * problem.F_map(Xp1,Xp2);
                    
                    %Diffusion matrix
                    K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) =...
                        K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain))...
                        * G' * B(end-problem.p:end) * wGPXIGA(iGP) * inverseJacobianX_Xi * problem.F_map(Xp1,Xp2);
                    
                    %% Integrate XIGA block
                    %extrnal heat source
                    f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        problem.rhs(globalGP, time) * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1,Xp2);
                    
                    %Capacity matrix
                    M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain))...
                        * ( F' * F ) * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1,Xp2);
                    
                    %Diffusion matrix
                    K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain))...
                        * ( G' * G ) * wGPXIGA(iGP) * inverseJacobianX_Xi * problem.F_map(Xp1,Xp2);
                    break;
                end
            end
            
        end
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
    else
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            localCoords = mapParentToLocal(rGP(iGP), Xp1, Xp2);
            globalCoords = mapParentToGlobal(rGP(iGP), Xp1, Xp2, problem, e);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);
            
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            inverseJacobianX_Xi = 1 / detJacobianX_Xi;
            
            %% Integrate IGA block
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGP(iGP), evaluateTemperature(e, localCoords, problem, solutionCoefficients), ...
                evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients))...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            %Diffusion matrix
            K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, evaluateTemperature(e, localCoords, problem, solutionCoefficients))...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGP(iGP) * inverseJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
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
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(localCoords(k),problem.p, problem.knotVector);
end

for i=1:length(x)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end

projectedCoefficients = projectionOperator(:, problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));

end


function [ projectedCoefficients ] = evaluateTemperatureSubElements(integrationDomainIndex, x, problem, solutionCoefficients,...
    modes, subDomainShapeFunctionCoefficients, shapeFunctionCoefficients, indexLocalEnrichedNodes, element, integrationDomain)
% EVALUATETEMPERATURESUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, problem.p+1 + modes * length(indexLocalEnrichedNodes) );

% N = zeros(length(x), problem.p+1);
N = zeros(length(x), problem.p + 1);

F = zeros(length(x), modes * numel(indexLocalEnrichedNodes));

parentCoords = x;

XiParametric1 = problem.knotVector( element + problem.p );
XiParametric2 = problem.knotVector( element + 1 + problem.p );

Xi1 = integrationDomain( integrationDomainIndex );
Xi2 = integrationDomain( integrationDomainIndex + 1 );

for k=1:length(x)

    [N(k,:), ~] = shapeFunctionsAndDerivativesIGASubElements(parentCoords(k), integrationDomainIndex,...
        subDomainShapeFunctionCoefficients, XiParametric1, XiParametric2, Xi1, Xi2, problem);
    
%     [N(k,:), ~] = BsplinesShapeFunctionsAndDerivatives(mapParentToLocal(parentCoords(k),...
%         XiParametric1, XiParametric2), problem.p, problem.knotVector);
    
    [F(k,:), ~] = PODModesAndDerivativesIGA(problem, mapParentToLocal(parentCoords(k),...
        XiParametric1, XiParametric2), modes,...
        problem.reductionOperator, shapeFunctionCoefficients, integrationDomainIndex, indexLocalEnrichedNodes, element);
end

for i=1:length(parentCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
end


projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * solutionCoefficients(problem.N:end);

end