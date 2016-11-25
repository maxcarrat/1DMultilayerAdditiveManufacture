function [ Phi, PhiDerivative ] = PODModesAndDerivativesIGA( problem, parametricCoordinates, numberOfModes,...
    PODCoefficients, coefficients, integrationSubDomainIndex, indexLocalEnrichedNodes, e, knotVector )
% PODMODESANDDERIVATIVESIGA Evaluate the POD modes and their derivatives

%% Matrix of modal shape functions
% Evaluate the enriched basis function Phi(x) using the approach presented in Fries
% and Belytschko "The extended/generalized finite element method:
% An overview of the method and its applications"(2010).

% Define Phi
Phi =  [];
% Define PhiDerivative
PhiDerivative = [];

for nodalIndex = 1:length(indexLocalEnrichedNodes)
    for iMode=1:numberOfModes
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        Xp1 = knotVector(e + problem.p);
        Xp2 = knotVector(e + problem.p + 1);
        
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        
        parentCoords = mapParametricToParent( parametricCoordinates, Xp1, Xp2);
        localSubElementCoords = mapGlobalToLocal( parentCoords, Xi1, Xi2);
        
        %map backward to the [-1,+1] Gauss integration domian
        mapIntegrationSubDomainBackward = 2/(Xi2 - Xi1);
        
        N_subElement =  0.5 * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        B_subElement =  0.5 * [-1, 1] * mapIntegrationSubDomainBackward;        % element shape function
        
        [N_Iga, B_Iga] = BsplinesShapeFunctionsAndDerivatives(...
            parametricCoordinates, problem.p, knotVector);
        
        %map backward to the global domain
        mapIntegrationParametricDomainBackward = 1 / problem.coords(end) - problem.coords(1);
        mapIntegrationDomainBackward = problem.B_map(Xp1, Xp2);
                
        %POD coefficients vector
        if indexLocalEnrichedNodes(nodalIndex) == problem.IGAdof %length( PODCoefficients ) - problem.p
            Phi_Nodal = PODCoefficients(end, iMode);
        else
            Phi_Nodal = PODCoefficients(1, iMode);
        end
        
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode)  - Phi_Nodal;
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode)  - Phi_Nodal;
        Phi_coeff = [Phi1, Phi2];
        
        Phi_iMode =  N_subElement * Phi_coeff';
        Phi_iModeContinuousDer =  mapIntegrationDomainBackward * ...
            B_subElement * Phi_coeff' * mapIntegrationParametricDomainBackward;
        
        derivative_1 = N_Iga(indexLocalEnrichedNodes(nodalIndex)) * Phi_iModeContinuousDer;
        derivative_2 =  B_Iga(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode * mapIntegrationParametricDomainBackward;
        
        mode = N_Iga(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode;
        
        %Phi modal basis
        Phi = [Phi, mode ];
        
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, derivative_1 + derivative_2];
    end   %end modal loop
end %end nodal index loop
end