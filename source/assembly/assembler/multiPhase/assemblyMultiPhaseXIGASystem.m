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

% On active elements use the refined domain as integration domain

refinedControlPoints = length(problem.reductionOperator)-1;
integrationDomain = linspace(-1, 1, ceil(refinedControlPoints/problem.XN));

localRefCP = refinedControlPoints/problem.XN;

for e=1:problem.N
    
    ldof = problem.p + 1;
    
    Xp1 = problem.knotVector(e + problem.p);
    Xp2 = problem.knotVector(e + problem.p + 1);
    
    detJacobianX_Xi = problem.coords(end) - problem.coords(1);
    
    inverseJacobianX_Xi = 1 / detJacobianX_Xi;
    detJacobianParameterToGlobal = problem.F_map(Xp1, Xp2);
    detJacobianTotal = detJacobianX_Xi * detJacobianParameterToGlobal;
    
    if e > (problem.N - problem.XN)

        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            localCoords = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
            globalCoords = mapParametricToGlobal(localCoords, problem);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
            
            B = B * inverseJacobianX_Xi;
            
            %Evaluate temperature@GP
            temperatureAtGaussPoint = evaluateTemperature(e, localCoords, problem, solutionCoefficients);
            previousTemperatureAtGaussPoint = evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients);
            
            
            %% Integrate IGA block
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGPXIGA(iGP) * detJacobianTotal;
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGPXIGA(iGP), temperatureAtGaussPoint, previousTemperatureAtGaussPoint )...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * detJacobianTotal;
            
            %Diffusion matrix
            K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, temperatureAtGaussPoint)...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * detJacobianTotal;
            
            elementEnrichedIndex = e - (problem.N - problem.XN);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = problem.IGAdof;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            modalDofs = length(indexLocalEnrichedNodes)*modes;
            
            for integrationSubDomainIndex =...
                    1 : ceil(localRefCP)-1
                
                if rGPXIGA(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXIGA(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*floor(localRefCP)+1:(elementEnrichedIndex-1)*floor(localRefCP)...
                        + ceil(localRefCP),:);

                    [F, G] = PODModesAndDerivativesIGA( problem, localCoords, modes, PODCoefficients,...
                        integrationDomain, integrationSubDomainIndex, indexLocalEnrichedNodes,...
                        e, problem.knotVector );
                                        
                    %Evaluate temperature@GP
                    temperatureAtGaussPoint = evaluateTemperatureSubElements(...
                        integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        solutionCoefficients, problem.modes,...
                        integrationDomain, indexLocalEnrichedNodes, e);
                    
                    previousTemperatureAtGaussPoint = evaluateTemperatureSubElements(...
                        integrationSubDomainIndex, rGPXIGA(iGP), problem,...
                        oldSolutionCoefficients, problem.modes,...
                        integrationDomain, indexLocalEnrichedNodes, e);
                    
                    %% Integrate Coupling block
                    %Capacity matrix
                    M_Coupling(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LMC(elementEnrichedIndex, :)) +...
                        problem.heatCapacity(globalCoords, temperatureAtGaussPoint, previousTemperatureAtGaussPoint)...
                        * F' * N(problem.LM(e, 1:ldof)) * wGPXIGA(iGP) *  detJacobianTotal;
                    
                    %Diffusion matrix
                    K_Coupling(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LMC(elementEnrichedIndex, :)) =...
                        K_Coupling(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LMC(elementEnrichedIndex, :)) +...
                        problem.k(globalCoords, time, temperatureAtGaussPoint)...
                        * G' * B(problem.LM(e, 1:ldof)) * wGPXIGA(iGP) * detJacobianTotal;
                    
                    %% Integrate XIGA block
                    %extrnal heat source
                    f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        problem.rhs(globalCoords, time) * wGPXIGA(iGP) * detJacobianTotal;
                    
                    %Capacity matrix
                    M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.heatCapacity(globalCoords, temperatureAtGaussPoint, previousTemperatureAtGaussPoint)...
                        * ( F' * F ) * wGPXIGA(iGP) * detJacobianTotal;
                    
                    %Diffusion matrix
                    K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        K_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        problem.k(globalCoords, time, temperatureAtGaussPoint)...
                        * ( G' * G ) * wGPXIGA(iGP) * detJacobianTotal;
                    break;
                end
            end
            
        end        
    else
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            localCoords = mapParentToLocal(rGP(iGP), Xp1, Xp2);
            globalCoords = mapParametricToGlobal(localCoords, problem);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);

            B = B * inverseJacobianX_Xi;
            
            %Evaluate temperature@GP
            temperatureAtGaussPoint = evaluateTemperature(e, localCoords, problem, solutionCoefficients);
            previousTemperatureAtGaussPoint = evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients);
            
            %% Integrate IGA block
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
                time) * wGP(iGP) * detJacobianTotal;
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.heatCapacity( rGP(iGP), temperatureAtGaussPoint, previousTemperatureAtGaussPoint)...
                * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianTotal;
            
            %Diffusion matrix
            K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                problem.k(globalCoords, time, temperatureAtGaussPoint)...
                * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianTotal;
            
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

N = zeros(length(x), numel(problem.coords));

localCoords = x;

for k=1:length(x)
    projectionOperator(k,1:size(N,2)) = BsplinesShapeFunctionsAndDerivatives(localCoords(k),problem.p, problem.knotVector);
end

projectedCoefficients = projectionOperator(:, problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));

end


function [ projectedCoefficients ] = evaluateTemperatureSubElements(integrationDomainIndex, x, problem,...
    solutionCoefficients, modes, shapeFunctionCoefficients, indexLocalEnrichedNodes, element)
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

parentCoords = x;

XiParametric1 = problem.knotVector( element + problem.p );
XiParametric2 = problem.knotVector( element + 1 + problem.p );

for k=1:numberOfProjectionPoints
    
    parametricCoords = mapParentToLocal(parentCoords(k), XiParametric1, XiParametric2);
    
    [NIga, ~] = BsplinesShapeFunctionsAndDerivatives(parametricCoords, problem.p, problem.knotVector);
    
    projectionOperator(k,1:size(problem.LM,2)+modes * numel(indexLocalEnrichedNodes)) = [NIga(problem.LM(element,:)), PODModesAndDerivativesIGA(problem, parametricCoords, modes,...
        problem.reductionOperator, shapeFunctionCoefficients, integrationDomainIndex,...
        indexLocalEnrichedNodes, element, problem.knotVector)];
end

projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * solutionCoefficients(problem.N:end);

end