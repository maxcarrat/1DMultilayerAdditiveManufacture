function [ Phi, PhiDerivative ] = PODModesAndDerivativesIGA( problem, localCoordinates, numberOfModes,...
    PODCoefficients, coefficients, integrationSubDomainIndex, indexLocalEnrichedNodes, e )
% PODMODESANDDERIVATIVESIGA Evaluate the POD modes and their derivatives

%% Matrix of modal shape functions
% Evaluate the enriched basis function Phi(x) using the approach presented in Fries
% and Belytschko "The extended/generalized finite element method:
% An overview of the method and its applications"(2010).

% Define Phi
Phi =  [];
for nodalIndex = 1:length(indexLocalEnrichedNodes)
    
    for iMode=1:numberOfModes
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        localSubElementCoords = mapGlobalToLocal( localCoordinates, Xi1, Xi2);
        
        N_subElement =  0.5 * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        
        % element shape function
        Xi1 = problem.knotVector( e + problem.p);
        Xi2 = problem.knotVector( e + problem.p + 1);
        [N, ~] = BsplinesShapeFunctionsAndDerivatives(...
            mapParentToLocal(localCoordinates, Xi1, Xi2), problem.p, problem.knotVector);
        
        %POD coefficients vector
        if indexLocalEnrichedNodes(nodalIndex) == 2
            Phi_Nodal = PODCoefficients(end, iMode);
        else
            Phi_Nodal = PODCoefficients(1, iMode);
        end
        
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode) - Phi_Nodal;
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode) - Phi_Nodal;
        Phi_coeff = [Phi1, Phi2];
        
        Phi_iMode = N_subElement * Phi_coeff';
        
        Phi_LocalSupport = Phi_iMode;
        
        %Phi modal basis
        Phi = [Phi, N(indexLocalEnrichedNodes(nodalIndex)) * Phi_LocalSupport ];
    end
    
end
%% Matrix of modes derivatives
% Evaluate the enriched basis function derivatives PhiDerivative(x)

% Define PhiDerivative
PhiDerivative = [];

for nodalIndex = 1:length(indexLocalEnrichedNodes)
    
    for iMode=1:numberOfModes
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        localSubElementCoords = mapGlobalToLocal( localCoordinates, Xi1, Xi2);
        
        mapIntegrationDomainForward =1; %(Xi2 - Xi1)/2 * 2 / problem.XN / 2;
        mapIntegrationDomainBackward = 2/(Xi2 - Xi1);
        
        N_subElement =   0.5...
            * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        B_subElement =  0.5 * [-1, 1];
        
        % element shape function and derivatives
        [N, B] = shapeFunctionsAndDerivatives(localCoordinates);
        
        %POD coefficients vector
        if indexLocalEnrichedNodes(nodalIndex) == 2
            Phi_Nodal = PODCoefficients(end, iMode);
        else
            Phi_Nodal = PODCoefficients(1, iMode);
        end
        
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode)  - Phi_Nodal;
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode)  - Phi_Nodal;
        Phi_coeff = [Phi1, Phi2];
        
%         Phi1Der = PODCoefficients(integrationSubDomainIndex, iMode);
%         Phi2Der = PODCoefficients(integrationSubDomainIndex+1, iMode);
%         Phi_coeff_der = [Phi1Der, Phi2Der];    
        
        Phi_iMode = mapIntegrationDomainForward * N_subElement * Phi_coeff';
        Phi_iModeContinuousDer = mapIntegrationDomainBackward * B_subElement * Phi_coeff';
        
        derivative_1 = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iModeContinuousDer;
        derivative_2 =  B(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode;
        
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, derivative_1 + derivative_2];
    end
end
end

