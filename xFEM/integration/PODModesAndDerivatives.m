function [ Phi, PhiDerivative ] = PODModesAndDerivatives( localCoordinates, numberOfModes,...
    PODCoefficients, coefficients, integrationSubDomainIndex, indexLocalEnrichedNodes )
% PODMODESANDDERIVATIVES Evaluate the POD modes and their derivatives
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

if length(indexLocalEnrichedNodes) == 1
    for iMode=1:numberOfModes
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        N_subElement =  0.5 * [(1 - localCoordinates)', (1 + localCoordinates)'];
        
        % coefficients of the element shape function
        N1 = coefficients(integrationSubDomainIndex);
        N2 = coefficients(integrationSubDomainIndex+1);
        N_coeff = [N1, N2];
        
        % element shape function

        N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];
        
        % POD coefficients vector
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode);
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode);
        Phi_coeff = [Phi1, Phi2];
        
        Phi_iMode = N_subElement * Phi_coeff';
        
        Phi_Nodal = PODCoefficients(end, iMode);
        
        Phi_LocalSupport = Phi_iMode - Phi_Nodal;
        
        % Phi modal basis
        Phi = [Phi, N(2) * Phi_LocalSupport];
    end
else
    for iMode=1:numberOfModes
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        N_subElement =  0.5 * [(1 - localCoordinates)', (1 + localCoordinates)'];
        
        %coefficients of the element shape function
        N1 = coefficients(integrationSubDomainIndex);
        N2 = coefficients(integrationSubDomainIndex+1);
        N_coeff = [N1, N2];
        
        %element shape function
        N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];
        
        %POD coefficients vector
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode);
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode);
        Phi_coeff = [Phi1, Phi2];
        
        Phi_iMode = N_subElement * Phi_coeff';
        
        Phi_Nodal = PODCoefficients(end, iMode);
        
        Phi_LocalSupport = abs(Phi_iMode - Phi_Nodal);
        
        %Phi modal basis
        Phi = [Phi, N * Phi_LocalSupport];
    end
    
end


%% Matrix of modes derivatives
% Evaluate the enriched basis function derivatives PhiDerivative(x)

% Define PhiDerivative
PhiDerivative = [];
if length(indexLocalEnrichedNodes) == 1
    for iMode=1:numberOfModes
        % Evaluate the ithMode basis function derivative PhiDerivative(x)
        % linear interpolation function on the integration subdomain
        B_subElement =  0.5 * [-1, 1];
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        N_subElement =  0.5 * [(1 - localCoordinates)', (1 + localCoordinates)'];
        
        %coefficients of the element shape function
        N1 = coefficients(integrationSubDomainIndex);
        N2 = coefficients(integrationSubDomainIndex+1);
        N_coeff = [N1, N2];
        
        %element shape function
        N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];
       
        %element shape function derivatives
        B = [B_subElement * ( 1 - N_coeff)', B_subElement * N_coeff'];
        
        %POD coefficients vector
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode);
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode);
        Phi_coeff = [Phi1, Phi2];
        
        Phi_Nodal = PODCoefficients(end, iMode);
        Phi_iMode = N_subElement * ( Phi_coeff - Phi_Nodal)';
        Phi_iModeDerivative = B_subElement * ( Phi_coeff - Phi_Nodal)';
   
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, B(2) * Phi_iMode + N(2) * Phi_iModeDerivative];
        
    end
else
    for iMode=1:numberOfModes
         % Evaluate the ithMode basis function derivative PhiDerivative(x)
        % linear interpolation function on the integration subdomain
        B_subElement =  0.5 * [-1, 1];
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        N_subElement =  0.5 * [(1 - localCoordinates)', (1 + localCoordinates)'];
        
        %coefficients of the element shape function
        N1 = coefficients(integrationSubDomainIndex);
        N2 = coefficients(integrationSubDomainIndex+1);
        N_coeff = [N1, N2];
        
        %element shape function
        N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];
       
        %element shape function derivatives
        B = [B_subElement * ( 1 - N_coeff)', B_subElement * N_coeff'];
        
        %POD coefficients vector
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode);
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode);
        Phi_coeff = [Phi1, Phi2];
        
        Phi_iMode = N_subElement * Phi_coeff';
        Phi_iModeDerivative = B_subElement * Phi_coeff';
   
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, B * Phi_iMode + N * Phi_iModeDerivative];
    end
end
end

