function [M, K, f] = assemblyXFEMGaussIntegrationSystem(problem, time, integrationOrder, integrationModalOrder, solutionCoefficients)
%   [M, K, f] = ASSEMBLY(problem) assembles the mass and the conductivity matrix and load vector
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

% On active elements use the refined domain as integration domain
refinedNodes = 2^problem.refinementDepth + 1;
integrationDomain = linspace(-1, 1, refinedNodes);
integrationCoefficients = linspace(-1, 1, refinedNodes);
subDomainShapeFunctionCoefficients = linspace(0, 1, refinedNodes);

for e=1:problem.N
    
    ldof = 2;
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    if e > (problem.N - problem.XN)
        
        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            [N, B] = shapeFunctionsAndDerivatives(rGPXFEM(iGP));
            
            %% Integrate FEM block
            %extrnal heat source
            f_FEM(problem.LM(e,1:ldof)) = f_FEM(problem.LM(e,1:ldof)) + N' * problem.rhs(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                time) * wGPXFEM(iGP) * problem.F_map(X1,X2);
            
            %Capacity matrix
            M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity * (N' * N) * wGPXFEM(iGP) * problem.F_map(X1,X2);
            
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
            
            for integrationSubDomainIndex = 1:numel(integrationDomain) - 1
                if rGPXFEM(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXFEM(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    Xi1 = integrationDomain(integrationSubDomainIndex);
                    Xi2 = integrationDomain(integrationSubDomainIndex + 1);
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*(refinedNodes - 1) + 1: (elementEnrichedIndex-1)*(refinedNodes - 1) + refinedNodes, :);
                    
                    mapIntegrationDomainForward = (Xi2 - Xi1)/2;
                    mapIntegrationDomainBackward = 2/(Xi2 - Xi1);
                    
                    %                     [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
                    [N, B] = shapeFunctionsAndDerivativesSubElements(rGPXFEM(iGP), integrationSubDomainIndex,...
                        subDomainShapeFunctionCoefficients);
                    
                    [F, G] = PODModesAndDerivativesGaussIntegration( rGPXFEM(iGP), modes, PODCoefficients,...
                        integrationCoefficients, integrationSubDomainIndex, indexLocalEnrichedNodes );
                    
                    globalGP = mapLocalToGlobal(rGPXFEM(iGP), X1, X2);
                    
                    %% Integrate Coupling block
                    %Capacity matrix
                    M_Coupling(((elementEnrichedIndex-1)*modalDofs + 1):modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(((elementEnrichedIndex-1)*modalDofs + 1):modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        problem.heatCapacity * F' * N * wGPXFEM(iGP) *  problem.F_map(X1,X2);
                    
                    %Diffusion matrix
                    K_Coupling(((elementEnrichedIndex-1)*modalDofs + 1):modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        K_Coupling(((elementEnrichedIndex-1)*modalDofs + 1):modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                        time, evaluateTemperatureSubElements(elementEnrichedIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * G' * B * wGPXFEM(iGP);
                    
                    %% Integrate XFEM block
                    %extrnal heat source
                    f_XFEM(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XFEM(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        problem.rhs(globalGP, time) * wGPXFEM(iGP) * problem.F_map(X1,X2);
                    
                    %Capacity matrix
                    M_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.heatCapacity * ( F' * F ) * wGPXFEM(iGP) * problem.F_map(X1,X2);
                    
                    %Diffusion matrix
                    K_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        K_XFEM(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGPXFEM(iGP), X1, X2),...
                        time, evaluateTemperatureSubElements(elementEnrichedIndex, rGPXFEM(iGP), problem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients, integrationCoefficients, indexLocalEnrichedNodes))...
                        * ( G' * G ) * wGPXFEM(iGP);
                    
                end
            end
            
            
        end
        
    else
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
            
            %% Integrate FEM block
            %extrnal heat source
            f_FEM(problem.LM(e,1:ldof)) = f_FEM(problem.LM(e,1:ldof)) + N' * problem.rhs(mapLocalToGlobal(rGP(iGP), X1, X2),...
                time) * wGP(iGP) * problem.F_map(X1,X2);
            
            %Capacity matrix
            M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity * (N' * N) * wGP(iGP) * problem.F_map(X1,X2);
            
            %Diffusion matrix
            K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_FEM(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
                time, evaluateTemperature(e, rGP(iGP), problem, solutionCoefficients)) * (B' * B) * wGP(iGP);
            
        end
    end
    
end

K = [K_FEM, K_Coupling(1:modalDofs,:)'; K_Coupling(1:modalDofs,:), K_XFEM(1:modalDofs,1:modalDofs)];
M = [M_FEM, M_Coupling(1:modalDofs,:)'; M_Coupling(1:modalDofs,:), M_XFEM(1:modalDofs,1:modalDofs)];
f = [f_FEM; f_XFEM(1:modalDofs)];


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
    modes, subDomainShapeFunctionCoefficients, shapeFunctionCoefficients,indexLocalEnrichedNodes)
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

F = zeros(length(x), modes);

localCoords = x;

for k=1:length(x)
    
    [N(k,:), ~] = shapeFunctionsAndDerivativesSubElements(localCoords(k), e,...
        subDomainShapeFunctionCoefficients);
    
    [F(k,:), ~] = PODModesAndDerivativesGaussIntegration(localCoords(k), modes,...
        problem.reductionOperator, shapeFunctionCoefficients, e, indexLocalEnrichedNodes);
end

for i=1:length(localCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
end


projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);

end