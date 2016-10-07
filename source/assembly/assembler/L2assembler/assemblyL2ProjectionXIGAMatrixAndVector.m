function [ M, f ] = assemblyL2ProjectionXIGAMatrixAndVector( solutionCoefficients, problem, integrationOrder,...
    integrationModalOrder, layerLength, initialTemperature, previousProblem )
%ASSEMBLYL2PROJECTIONMATRIXANDVECTOR assemble the matrix M and the vector f for
%the L2 projection system
%Input:
%solutionCoefficients = coefficients of the solutions
%problem = struct which defines the initial boundary value problem
%integrationOrder = number of quadrature points
%integrationModalOrder = numeber of quadratur points in enriched elements
%layerLength = length of bar layers
%initialTemperature = initial temperature of the system
%previousproblem = struct of the poisson problem in the previous layer
%Output:
%M = L2 projection matrix integral(N'N detJx detJxi)
%f = L2 rhs integral(N' functionToProject detJx detJxi)

%% Allocate matrices

%IGA  matrix
M_IGA = zeros(problem.IGAdof,problem.IGAdof);

%IGA  vector
f_IGA = zeros(problem.IGAdof, 1);

%XIGA  matrix
M_XIGA = zeros(problem.XIGAdof,problem.XIGAdof);

%XIGA  vector
f_XIGA = zeros(problem.XIGAdof, 1);

%Coupling  matrix
M_Coupling = zeros(problem.XIGAdof,problem.IGAdof);

%Gauss integration
[rGP, wGP] = gaussPoints( integrationOrder );
[rGPXIGA, wGPXIGA] = gaussPoints( integrationModalOrder );

numberOfIntegrationPoints = length(rGP);
numberOfModalIntegrationPoints = length(rGPXIGA);

modes = problem.modes;
incrementModes = 0;

%On active elements use the refined domain as integration domain

refinedControlPoints = 2^problem.refinementDepth + problem.p;
integrationDomain = linspace(-1, 1, ceil(refinedControlPoints/problem.XN));
subDomainShapeFunctionCoefficients = linspace(0, 1, ceil(refinedControlPoints/problem.XN));

%loop over elements
for e=1:problem.N
    
    ldof = problem.p + 1;
    
    %left right end of the knot span
    Xp1 = problem.knotVector(e + problem.p);
    Xp2 = problem.knotVector(e + problem.p + 1);
    
    %% Enriched element
    if e > problem.N - problem.XN    %if element of the new mesh is enriched
        
        %In this element the temperature values are constantly equal to the
        %initial temperature, but they are projected onto an enriched
        %element.
        
        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            %map GP onto the parametric space
            localCoords = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
            
            %map GP onto the global space
            globalCoords = mapParentToGlobal(rGPXIGA(iGP), Xp1, Xp2, problem, e);
            
            %evaluate BSplines at integration points
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
            temperatureAtGaussPoint = evaluateTemperature(e, globalCoords, localCoords, problem,...
                solutionCoefficients, layerLength, initialTemperature);
            
            %Jacobian from parametric to global and from parent to
            %parametric
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            detJacobianParameterToGlobal = problem.F_map(Xp1, Xp2);
            
            %% Integrate IGA block
            
            % vector
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                * wGPXIGA(iGP) * detJacobianX_Xi * detJacobianParameterToGlobal;
            
            % matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * detJacobianX_Xi * detJacobianParameterToGlobal;
            
            elementEnrichedIndex = e - (problem.N - problem.XN);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = problem.IGAdof;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            modalDofs = length(indexLocalEnrichedNodes)*modes;
            
            %loop over sub-domains
            for integrationSubDomainIndex =...
                    1 : ceil(refinedControlPoints/problem.XN)-1
                
                %if the GP is inside the sub-domain
                if rGPXIGA(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXIGA(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    PODCoefficients = problem.reductionOperator(...
                        (elementEnrichedIndex-1)*floor(refinedControlPoints/problem.XN)+1:(elementEnrichedIndex-1)*floor(refinedControlPoints/problem.XN)...
                        + ceil(refinedControlPoints/problem.XN),:);
                    
                    %map GP onto parametric space
                    localCoords = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
                    
                    %map GP onto global space
                    globalCoords = mapParentToGlobal(rGPXIGA(iGP), Xp1, Xp2, problem, e);
                    
                    %evaluate BSplines and POD-modal enrichment functions
                    %@GP
                    [N, ~] = shapeFunctionsAndDerivativesIGASubElements( rGPXIGA(iGP), integrationSubDomainIndex,...
                        subDomainShapeFunctionCoefficients,  Xp1, Xp2,...
                        integrationDomain(integrationSubDomainIndex), integrationDomain(integrationSubDomainIndex+1), problem );
                    [F, ~] = PODModesAndDerivativesIGA( problem, localCoords, modes, PODCoefficients,...
                        integrationDomain, integrationSubDomainIndex, indexLocalEnrichedNodes, e, problem.knotVector );
                    
                    %Evaluate temperature@GP
                    temperatureAtGaussPoint = evaluateTemperature(e, globalCoords, localCoords, problem,...
                        solutionCoefficients, layerLength, initialTemperature);
                    
                    
                    %% Integrate Coupling block
                    % matrix
                    M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) =...
                        M_Coupling(((elementEnrichedIndex-2)*incrementModes + 1):...
                        ((elementEnrichedIndex-2)*incrementModes) + modalDofs, problem.LMC(elementEnrichedIndex, :)) +...
                        F' * N(end-problem.p:end) * wGPXIGA(iGP) *  detJacobianX_Xi * detJacobianParameterToGlobal;
                    
                    %% Integrate XIGA block
                    % vector
                    f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) = f_XIGA(problem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                        temperatureAtGaussPoint * wGPXIGA(iGP) * detJacobianX_Xi * detJacobianParameterToGlobal;
                    
                    % matrix
                    M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                        M_XIGA(problem.LME(elementEnrichedIndex, 1:modalDofs), problem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                        ( F' * F ) * wGPXIGA(iGP) * detJacobianX_Xi * detJacobianParameterToGlobal;
                    
                    break;
                end
            end
            
        end
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
        %% Previously enriched element
    elseif e >= (problem.N - problem.XN) && e <= (problem.N - problem.XN)
        
        %In this element the temperature values, i.e. the function
        %values to be projected, are evaluated on the enriched element
        %of the previous mesh. It means the shape functions are the
        %normal BSplines of the non-enriched Iga elements but the
        %temperature is evaluated onto an enriched element.
        
        
        % Gauss integration
        for iGP = 1:numberOfModalIntegrationPoints
            
            elementEnrichedIndex = e - (problem.N - problem.XN - 1);
            
            if elementEnrichedIndex == 1
                indexLocalEnrichedNodes = previousProblem.IGAdof;
            else
                indexLocalEnrichedNodes = [1, 2];
            end
            
            
            for integrationSubDomainIndex =...
                    1 : ceil(refinedControlPoints/problem.XN)-1
                
                if rGPXIGA(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                        rGPXIGA(iGP) > integrationDomain(integrationSubDomainIndex)
                    
                    %map GP onto the global space
                    globalCoords = mapParentToGlobal(rGPXIGA(iGP), Xp1, Xp2, problem, e);
                    
                    %map GP onto the parametric space
                    localCoords = mapParentToLocal(rGPXIGA(iGP), Xp1, Xp2);
                    
                    %evaluate shape function at GP and Jacobian
                    [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);
                    JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
                    detJacobianX_Xi = norm(JacobianX_Xi);
                    
                    %Evaluate temperature@GP
                    temperatureAtGaussPoint = evaluateTemperatureSubElements(integrationSubDomainIndex,...
                        rGPXIGA(iGP), globalCoords, problem, previousProblem,...
                        solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients,...
                        integrationDomain, indexLocalEnrichedNodes, e, integrationDomain, layerLength, initialTemperature);
                    
                    
                    %% Integrate IGA block
                    % vector
                    f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                        * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
                    
                    %Capacity matrix
                    M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                        (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGPXIGA(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
                    
                    break;
                end
            end
            
        end
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
        %% Non-enriched element
    else
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            localCoords = mapParentToLocal(rGP(iGP), Xp1, Xp2);
            globalCoords = mapParentToGlobal(rGP(iGP), Xp1, Xp2, problem, e);
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);
            
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            detJacobianX_Xi = norm(JacobianX_Xi);
            
            %Evaluate temperature@GP:
            temperatureAtGaussPoint = evaluateTemperature(e, globalCoords, localCoords, problem,...
                solutionCoefficients, layerLength, initialTemperature);
            
            
            %% Integrate IGA block
            
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xp1, Xp2);
            
        end
    end
    
end

%% Assembly
M = [M_IGA, M_Coupling'; M_Coupling, M_XIGA];
f = [f_IGA; f_XIGA];


end


function [ projectedCoefficients ] = evaluateTemperature(e, xGlobal, xLocal,...
    problem, solutionCoefficients, layerLength, initialTemperature)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(xLocal);

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(xLocal), 2);

localCoords = xLocal;

for k=1:length(xLocal)
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(localCoords(k),problem.p, problem.knotVector);
end

for i=1:length(xLocal)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end

if xGlobal < problem.coords(end) && xGlobal > ( problem.coords(end) - layerLength)
    projectedCoefficients = initialTemperature;
else
    projectedCoefficients = projectionOperator(:, problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));
end

end


function [ projectedCoefficients ] = evaluateTemperatureSubElements(integrationDomainIndex,...
    x, globalCoords, problem, previousProblem, solutionCoefficients,...
    modes, subDomainShapeFunctionCoefficients, shapeFunctionCoefficients,...
    indexLocalEnrichedNodes, element, integrationDomain, layerLength, initialTemperature)
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

N = zeros(length(x), problem.p + 1);

F = zeros(length(x), modes * numel(indexLocalEnrichedNodes));

parentCoords = x;

XiParametric1 = problem.knotVector( element + problem.p );
XiParametric2 = problem.knotVector( element + 1 + problem.p );

Xi1 = integrationDomain( integrationDomainIndex );
Xi2 = integrationDomain( integrationDomainIndex + 1 );

for k=1:length(x)
    
    [N(k,:), ~] = shapeFunctionsAndDerivativesIGASubElements(parentCoords(k), integrationDomainIndex,...
        subDomainShapeFunctionCoefficients, XiParametric1, XiParametric2, Xi1, Xi2, previousProblem);
    
    [F(k,:), ~] = PODModesAndDerivativesIGA(previousProblem, mapParentToLocal(parentCoords(k),...
        XiParametric1, XiParametric2), modes,...
        previousProblem.reductionOperator, shapeFunctionCoefficients, integrationDomainIndex,...
        indexLocalEnrichedNodes, element, previousProblem.knotVector);
end

for i=1:length(parentCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,1:end), F(i,:)];
end


if globalCoords <= problem.coords(end) && globalCoords > ( problem.coords(end) - layerLength)
    projectedCoefficients = initialTemperature;
else
    projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * solutionCoefficients(problem.N-1:end);
end


end