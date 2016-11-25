function [ Phi, PhiDerivative ] = PODModesAndDerivativesMultiscale( localCoordinates, elementGlobalCoords, numberOfModes,...
    PODCoefficients, coefficients, integrationSubDomainIndex, indexLocalEnrichedNodes )
%PODMODESANDDERIVATIVESMULTISCALE Evaluate the POD modes and their derivatives
%   localCoordinate = integration point local coordinate
%   numberOfModes = number of modal functions of the enriched problem
%   PODCoefficients = nodal values of POD basis
%   integrationSubDomainIndex = index of the integration sub-domain

%% Matrix of modal shape functions
% Evaluate the enriched basis function Phi(x) using the approach presented in Fries
% and Belytschko "The extended/generalized finite element method:
% An overview of the method and its applications"(2010).

% Define Phi
Phi =  [];
% Define PhiDerivative
PhiDerivative = [];

numberOfLocalEnrichedNodes = length(indexLocalEnrichedNodes);

for nodalIndex = 1:numberOfLocalEnrichedNodes
    
    for iMode=1:numberOfModes
        
        %Evaluate the ithMode basis function Phi(x)
        %linear interpolation function on the integration subdomain
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        localSubElementCoords = mapGlobalToLocal( localCoordinates, Xi1, Xi2);
        mapIntegrationSubDomainBackward = 2/(Xi2 - Xi1);
        mapIntegrationDomainBackward= 2/(elementGlobalCoords(2) - elementGlobalCoords(1));

        N_subElement =   0.5...
            * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        B_subElement =  0.5 * [-1, 1];
        
        
        % element shape function and derivatives
        [N, B] = shapeFunctionsAndDerivatives(localCoordinates);
        
        %POD coefficients vector
        if indexLocalEnrichedNodes(nodalIndex) == 2
            Phi_Nodal = PODCoefficients(end, iMode);
            
            Phi1 = PODCoefficients(integrationSubDomainIndex, iMode) - Phi_Nodal;
            Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode) - Phi_Nodal;
            Phi_coeff = [Phi1, Phi2];
            
            Phi_iMode = N_subElement * Phi_coeff';
            Phi_iModeContinuousDer = mapIntegrationSubDomainBackward ...
                * B_subElement * Phi_coeff' * mapIntegrationDomainBackward;
            
%             Phi_LocalSupport =  Phi_iMode;
            Phi_LocalSupport = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode;
            
            derivative_1 = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iModeContinuousDer;
            derivative_2 = B(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode *...
                mapIntegrationDomainBackward;            
            
        else
            Phi_Nodal = PODCoefficients(1, iMode);
            
            Phi1 = PODCoefficients(integrationSubDomainIndex, iMode) - Phi_Nodal;
            Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode) - Phi_Nodal;
            Phi_coeff = [Phi1, Phi2];
            
            Phi_iMode = N_subElement * Phi_coeff';
            Phi_iModeContinuousDer = mapIntegrationSubDomainBackward ...
                * B_subElement * Phi_coeff' * mapIntegrationDomainBackward;
            
            Phi_LocalSupport = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode;
            
            derivative_1 = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iModeContinuousDer;
            derivative_2 = B(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode *...
                mapIntegrationDomainBackward;
            
        end
        
        %Phi modal basis
        Phi = [Phi, Phi_LocalSupport ];
        
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, derivative_1 + derivative_2];
    end
    
end