function [ K_Coupling, K_XFEM, M_Coupling, M_XFEM, f_XFEM, modalDofs] = integrateExtendedTerms( problem, integrationOrder )
%INTEGRATEEXTENDEDTERMS pre integrate the extended terms for linear
% problems, assumig zero rhs

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


for e=1:problem.XN
    
    elementOffset = problem.N - problem.XN;
    X1 = problem.coords(e + elementOffset);
    X2 = problem.coords(e+1 + elementOffset);
    
    if e == 1 
        indexLocalEnrichedNodes = 2; %rhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    modes = problem.modes;
    modalDofs = length(indexLocalEnrichedNodes) * modes;
    
    % On active elements use the refined domain as integration domain
    refinedNodes = 2^problem.refinementDepth + 1;
    integrationDomain = linspace(-1, 1, refinedNodes);
    subDomainShapeFunctionCoefficients = linspace(0, 1, refinedNodes);
    
    for integrationSubDomain=1:length(integrationDomain)-1
        
        Xi1 = integrationDomain(integrationSubDomain);
        Xi2 = integrationDomain(integrationSubDomain + 1);
        
        PODCoefficients = problem.reductionOperator(...
            (e-1)*(refinedNodes - 1) + 1: (e-1)*(refinedNodes - 1) + refinedNodes, :);
        
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
                problem.B_map(X1,X2) * problem.k * G' * B * wGP(iGP) * mapIntegrationDomainBackward;
            
            %% Integrate XFEM block
            %extrnal heat source
            f_XFEM(problem.LME(e,1:modalDofs)) = f_XFEM(problem.LME(e,1:modalDofs)) +  F' *...
                0.0 * wGP(iGP) * problem.F_map(X1,X2) * mapIntegrationDomainForward;
            
            %Capacity matrix
            M_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) = M_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) +...
                problem.heatCapacity * ( F' * F ) * wGP(iGP) * problem.F_map(X1,X2) * mapIntegrationDomainForward;
            
            %Diffusion matrix
            K_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) = K_XFEM(problem.LME(e, 1:modalDofs), problem.LME(e, 1:modalDofs)) +...
                problem.B_map(X1,X2) * problem.k * ( G' * G ) * wGP(iGP) * mapIntegrationDomainBackward;
            
        end
    end
    
end


end

