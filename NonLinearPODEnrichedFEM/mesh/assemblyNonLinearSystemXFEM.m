function [M, K, f] = assemblyNonLinearSystemXFEM(problem, time, integrationOrder, solutionCoefficients)
%   ASSEMBLYNONLINEARSYSTEMXFEM assembles the mass and the conductivity matrix and load vector
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

numberOfIntegrationPoints = length(rGP);

modes = problem.modes;

for e=1:problem.N
    
    ldof = 2;
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
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
            time, evaluateTemperature(e, rGP(iGP), problem, solutionCoefficients))...
            * (B' * B) * wGP(iGP);
        
    end
    
end


for e=1:problem.XN
    
    elementOffset = problem.N - problem.XN;
    X1 = problem.coords(e + elementOffset);
    X2 = problem.coords(e+1 + elementOffset);
    
    if e == 1
        indexLocalEnrichedNodes = 2; %rhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    modalDofs = length(indexLocalEnrichedNodes)*modes;
    
    % On active elements use the refined domain as integration domain
    refinemedNodes = 2^problem.refinementDepth + 1;
    integrationDomain = linspace(-1, 1, refinemedNodes);
    subDomainShapeFunctionCoefficients = linspace(0, 1, refinemedNodes);
    
    for integrationSubDomain=1:length(integrationDomain)-1
        
        Xi1 = integrationDomain(integrationSubDomain);
        Xi2 = integrationDomain(integrationSubDomain + 1);
        
        PODCoefficients = problem.reductionOperator(...
            (e-1)*(refinemedNodes - 1) + 1: (e-1)*(refinemedNodes - 1) + refinemedNodes, :);
        
        mapIntegrationDomainForward = (Xi2 - Xi1)/2;
        mapIntegrationDomainBackward = 2/(Xi2 - Xi1);
        
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            [N, B] = shapeFunctionsAndDerivativesSubElements(rGP(iGP), integrationSubDomain,...
                subDomainShapeFunctionCoefficients);
            %                 [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
            
            [F, G] = PODModesAndDerivatives(rGP(iGP), modes, PODCoefficients,...
                subDomainShapeFunctionCoefficients, integrationSubDomain, indexLocalEnrichedNodes);
            
            globalGP = mapLocalToGlobal(mapLocalToGlobal(rGP(iGP), Xi1, Xi2), X1, X2);
            
            %% Integrate Coupling block
            %Capacity matrix
            M_Coupling(((e-1)*modalDofs + 1):modalDofs, problem.LMC(e, :)) = M_Coupling(((e-1)*modalDofs + 1):modalDofs, problem.LMC(e, :)) +...
                problem.heatCapacity * F' * N * wGP(iGP) *  problem.F_map(X1,X2) * mapIntegrationDomainForward;
            
            %Diffusion matrix
            K_Coupling(((e-1)*modalDofs + 1):modalDofs, problem.LMC(e, :)) = K_Coupling(((e-1)*modalDofs + 1):modalDofs, problem.LMC(e, :)) +...
                problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
                time, evaluateTemperatureSubElements(e, rGP(iGP), problem,...
                solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients)) * G' * B * wGP(iGP) * mapIntegrationDomainBackward;
            
            %% Integrate XFEM block
            %extrnal heat source
            f_XFEM(problem.LME(e,1:modalDofs)) = f_XFEM(problem.LME(e,1:modalDofs)) +  F' *...
                problem.rhs(globalGP, time) * wGP(iGP) * problem.F_map(X1,X2) * mapIntegrationDomainForward;
            
            %Capacity matrix
            M_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) = M_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) +...
                problem.heatCapacity * ( F' * F ) * wGP(iGP) * problem.F_map(X1,X2) * mapIntegrationDomainForward;
            
            %Diffusion matrix
            K_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) = K_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) +...
                problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
                time, evaluateTemperatureSubElements(e, rGP(iGP), problem,...
                solutionCoefficients, problem.modes, subDomainShapeFunctionCoefficients))...
                * ( G' * G ) * wGP(iGP) * mapIntegrationDomainBackward;
            
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

X1 = problem.coords(e);
X2 = problem.coords(e+1);


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
    modes, subDomainShapeFunctionCoefficients)
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
B = zeros(length(x), 2);

F = zeros(length(x), modes);
G = zeros(length(x), modes);

localCoords = x;

for k=1:length(x)
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivativesSubElements(localCoords(k), e,...
        subDomainShapeFunctionCoefficients);
    [F(k,:), G(k,:)] = PODModesAndDerivatives(localCoords(k), modes, problem.reductionOperator,...
        subDomainShapeFunctionCoefficients, e, 2);
end

for i=1:length(localCoords)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
end


projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);

end