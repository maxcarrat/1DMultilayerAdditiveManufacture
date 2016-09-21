function [M, K, f] = assemblyMultiPhaseXFEMSystem(problem, time, integrationOrder, integrationModalOrder,...
    solutionCoefficients, oldSolutionCoefficients)
%ASSEMBLYMULTIPHASEXFEMSYSTEM: assembles the mass and the conductivity matrix and load vector
%   problem = definition of the boundary value problem
%   time = current time
%   integrationOrder = number of integration points

%% Allocate matrices
% FEM block
%FEM conductivity matrix
K_FEM = zeros(problem.FEMdof,problem.FEMdof);
%FEM capacity matrix
M_FEM = zeros(problem.FEMdof,problem.FEMdof);
%FEM load vector
f_FEM = zeros(problem.FEMdof, 1);

% XFEM block
%XFEM conductivity matrix
K_XFEM = zeros(problem.XFEMdof,problem.XFEMdof);
%XFEM capacity matrix
M_XFEM = zeros(problem.XFEMdof,problem.XFEMdof);
%XFEM load vector
f_XFEM = zeros(problem.XFEMdof, 1);

% Coupling block
%Coupling conductivity matrix
K_Coupling = zeros(problem.XFEMdof,problem.FEMdof);
%Coupling capacity matrix
M_Coupling = zeros(problem.XFEMdof,problem.FEMdof);


%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
[rGPXFEM, wGPXFEM] = gaussPoints( integrationModalOrder );

numberOfIntegrationPoints = length(rGP);
numberOfModalIntegrationPoints = length(rGPXFEM);

modes = problem.modes;
incrementModes = 0;

% On active elements use the refined domain as integration domain
refinedNodes = 2^problem.refinementDepth + 1;
integrationDomain = linspace(-1, 1, ceil(refinedNodes/problem.XN));
integrationCoefficients = linspace(-1, 1, ceil(refinedNodes/problem.XN));
subDomainShapeFunctionCoefficients = linspace(0, 1, ceil(refinedNodes/problem.XN));

for e=1:problem.N
    
    ldof = 2;
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    if e > (problem.N - problem.XN)
        
        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            [N, B] = shapeFunctionsAndDerivatives(rGPXFEM(iGP));
            
            %% Integrate FEM block
            %External heat source
            f_FEM(problem.LM(e,1:ldof)) = f_FEM(problem.LM(e,1:ldof)) + N' * problem.rhs(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                time) * wGPXFEM(iGP) * problem.F_map(X1,X2);
            
            %Capacity matrix
            M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                evaluateTemperature(e, rGPXFEM(iGP), problem, solutionCoefficients),...
                evaluateTemperature(e, rGPXFEM(iGP), problem, oldSolutionCoefficients))...
                * (N' * N) * wGPXFEM(iGP) * problem.F_map(X1,X2);
            
            %Diffusion matrix
            K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                time, evaluateTemperature(e, rGPXFEM(iGP), problem, solutionCoefficients)) * (B' * B) * wGPXFEM(iGP);
            
            
            elementEnrichedIndex = e - (problem.N - problem.XN);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = 2;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            modalDofs = length(indexLocalEnrichedNodes)*modes;

            for integrationSubDomainIndex =...
                    1 : ceil(refinedNodes/problem.XN)-1
                    
                if rGPXFEM(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXFEM(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*floor(refinedNodes/problem.XN)+1:(elementEnrichedIndex-1)*floor(refinedNodes/problem.XN)...
                        + ceil(refinedNodes/problem.XN),:);

                    [N, B] = shapeFunctionsAndDerivativesSubElements(rGPXFEM(iGP), integrationSubDomainIndex,...
                        subDomainShapeFunctionCoefficients, integrationDomain);
                    
                    [F, G] = PODModesAndDerivativesGaussIntegration(rGPXFEM(iGP), modes, PODCoefficients,...
                        integrationCoefficients, integrationSubDomainIndex, indexLocalEnrichedNodes );
                    
                    globalGP = mapLocalToGlobal(rGPXFEM(iGP), X1, X2);
                    
                    %% Integrate Coupling block
                    %Capacity matrix
                    M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                    ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * F' * N * wGPXFEM(iGP) *  problem.F_map(X1,X2);
                    
                    %Diffusion matrix
                    K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) =...
                        K_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):((elementEnrichedIndex-2)*incrementModes) + modalDofs...
                        , problem.LMC(elementEnrichedIndex, :)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * G' * B * wGPXFEM(iGP) * problem.B_map(X1,X2);
                    
                    %% Integrate XFEM block
                    %extrnal heat source
                    f_XFEM(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XFEM(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        problem.rhs(globalGP, time) * wGPXFEM(iGP) * problem.F_map(X1,X2);
                    
                    %Capacity matrix
                    M_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.heatCapacity(globalGP,...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes),...
                        evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        oldSolutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * ( F' * F ) * wGPXFEM(iGP) * problem.F_map(X1,X2);
                    
                    %Diffusion matrix
                    K_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        K_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.k(globalGP,...
                        time, evaluateTemperatureSubElements(integrationSubDomainIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * ( G' * G ) * wGPXFEM(iGP) * problem.B_map(X1,X2);
                    break;
                end
            end
            
        end
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
    else
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            globalGP = mapLocalToGlobal(rGPXFEM(iGP), X1, X2);
            
            [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
            
            %% Integrate FEM block
            %extrnal heat source
            f_FEM(problem.LM(e,1:ldof)) = f_FEM(problem.LM(e,1:ldof)) + N' * problem.rhs(globalGP,...
                time) * wGP(iGP) * problem.F_map(X1,X2);
            
            %Capacity matrix
            M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity(globalGP,...
                evaluateTemperature(e, rGP(iGP), problem, solutionCoefficients),...
                evaluateTemperature(e, rGP(iGP), problem, oldSolutionCoefficients))...
                * (N' * N) * wGP(iGP) * problem.F_map(X1,X2);
            
            %Diffusion matrix
            K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.B_map(X1,X2) * problem.k(globalGP,...
                time, evaluateTemperature(e, rGP(iGP), problem, solutionCoefficients)) * (B' * B) * wGP(iGP);
            
        end
    end
    
end

K = [K_FEM, K_Coupling'; K_Coupling, K_XFEM];
M = [M_FEM, M_Coupling'; M_Coupling, M_XFEM];
f = [f_FEM; f_XFEM];


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
    modes, subDomainShapeFunctionCoefficients, shapeFunctionCoefficients, indexLocalEnrichedNodes)
% EVALUATETEMPERATURESUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes

numberOfProjectionPoints = length(x);
refinedNodes = 2^problem.refinementDepth + 1;
projectionOperator = zeros(numberOfProjectionPoints, length(solutionCoefficients(problem.N:end)) );
integrationDomain = linspace(-1, 1, ceil(refinedNodes/problem.XN));

N = zeros(length(x), 2);

F = zeros(length(x), modes * numel(indexLocalEnrichedNodes));

localCoords = x;

for k=1:length(x)
    
    [N(k,:), ~] = shapeFunctionsAndDerivativesSubElements(localCoords(k), e,...
        subDomainShapeFunctionCoefficients, integrationDomain);
    
    [F(k,:), ~] = PODModesAndDerivativesGaussIntegration(localCoords(k), modes,...
        problem.reductionOperator, shapeFunctionCoefficients, e, indexLocalEnrichedNodes);
end

for i=1:length(localCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
end


projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);

end