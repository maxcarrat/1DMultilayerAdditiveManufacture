function [ Phi, PhiDerivative ] = PODModesAndDerivativesIGA( problem, parametricCoordinates, numberOfModes,...
    PODCoefficients, coefficients, integrationSubDomainIndex, indexLocalEnrichedNodes, e )
% PODMODESANDDERIVATIVESIGA Evaluate the POD modes and their derivatives
%Input: 
%problem = poisson problem struct
%parametricCoordinates = parametric coordinates of the quadrature point
%numberOfModes = number of POD-modes
%PODCoefficients = coefficients of POD eigenvectors
%integrationSubDomainIndex = index of the actual integration sub-doamain
%indexLocalEnrichedNodes = node index 
%e = element index
%Output:
%Phi = eXtended basis functions
%PhiDerivative = eXtended basis function derivatives

%% Matrix of modal shape functions
% Evaluate the enriched basis function Phi(x) using the approach presented in Fries
% and Belytschko "The extended/generalized finite element method:
% An overview of the method and its applications"(2010).

% Define Phi
Phi =  [];

%loop over the enriched degrees of freedom
for nodalIndex = 1:length(indexLocalEnrichedNodes)
    
%loop over the POD-modes   
    for iMode=1:numberOfModes
        
        %% Evaluate the ithMode basis function Phi(x)
        %Linear interpolation function on the integration subdomain.
        %The modal basis are defined using linear polynomials defined over
        %the sub-domains, the sub-domains node coefficients are the
        %POD eigenvector coefficients.
        
        %knot span left ad right end
        Xp1 = problem.knotVector(e + problem.p);
        Xp2 = problem.knotVector(e + problem.p + 1);
        
        %integration sub-domain left and right end
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        
        %map from parametric to parent space
        parentCoords = mapParametricToParent( parametricCoordinates, Xp1, Xp2);
        
        %map from parent to sub-domain space where the linear interpolation
        %functions are defined
        localSubElementCoords = mapGlobalToLocal( parentCoords, Xi1, Xi2);
        
        %sub-domain linear interpolation functions
        N_subElement =  0.5 * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        
        % element BSplines
        [N, ~] = BsplinesShapeFunctionsAndDerivatives(...
            parametricCoordinates, problem.p, problem.knotVector);  
        
        %POD coefficients vector
        %taken from: Belytschko T, Moes N, Usui S, Parimi C. Arbitrary discontinuities in finite elements.
        %International Journal for Numerical Methods in Engineering 2001; 50:993–1013.
        if indexLocalEnrichedNodes(nodalIndex) == problem.IGAdof
            Phi_Nodal = PODCoefficients(end, iMode);
        else
            Phi_Nodal = PODCoefficients(1, iMode);
        end
        
        %generate eXtended basis coefficients
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode) - Phi_Nodal;
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode) - Phi_Nodal;
        Phi_coeff = [Phi1, Phi2];
        
        %apply partition of unity method
        Phi_iMode = N_subElement * Phi_coeff';        
        
        %Phi modal basis
        Phi = [Phi, N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode ];
    end
    
end
%% Matrix of modes derivatives
% Evaluate the enriched basis function derivatives PhiDerivative(x)

% Define PhiDerivative
PhiDerivative = [];

%loop over the enriched degrees of freedom
for nodalIndex = 1:length(indexLocalEnrichedNodes)
    
%loop over the POD-modes   
    for iMode=1:numberOfModes
        
        % Evaluate the ithMode basis function Phi(x)
        % linear interpolation function on the integration subdomain
        
        %knot span left ad right end
        Xp1 = problem.knotVector(e + problem.p);
        Xp2 = problem.knotVector(e + problem.p + 1);
        
        %integration sub-domain left and right end
        Xi1 = coefficients(integrationSubDomainIndex);
        Xi2 = coefficients(integrationSubDomainIndex + 1);
        
        %map from parametric to parent space
        parentCoords = mapParametricToParent( parametricCoordinates, Xp1, Xp2);
        
        %map from parent to sub-domain space where the linear interpolation
        %functions are defined
        localSubElementCoords = mapGlobalToLocal( parentCoords, Xi1, Xi2);
        
        %inverse of the Jacobian from parent to parametric
        mapIntegrationDomainBackward = 2/(Xp2 - Xp1);
        %inverse of the Jacobian from sub-domain to parent
        mapIntegrationSubDomainBackward = 2/(Xi2 - Xi1);
        
        %sub-domain linear interpolation functions and derivatives       
        N_subElement =   0.5...
            * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
        B_subElement =  0.5 * [-1, 1];
        
        % element shape function and derivatives
        [N, B] = BsplinesShapeFunctionsAndDerivatives(...
            parametricCoordinates, problem.p, problem.knotVector);
        
        %POD coefficients vector
        %taken from: Belytschko T, Moes N, Usui S, Parimi C. Arbitrary discontinuities in finite elements.
        %International Journal for Numerical Methods in Engineering 2001; 50:993–1013.
        if indexLocalEnrichedNodes(nodalIndex) == problem.IGAdof 
            Phi_Nodal = PODCoefficients(end, iMode);
        else
            Phi_Nodal = PODCoefficients(1, iMode);
        end
        
        %generate eXtended basis coefficients        
        Phi1 = PODCoefficients(integrationSubDomainIndex, iMode)  - Phi_Nodal;
        Phi2 = PODCoefficients(integrationSubDomainIndex+1, iMode)  - Phi_Nodal;
        Phi_coeff = [Phi1, Phi2];
        
        %apply partition of unity method        
        Phi_iMode =  N_subElement * Phi_coeff';
        Phi_iModeContinuousDer =  mapIntegrationDomainBackward * mapIntegrationSubDomainBackward * ...
             B_subElement * Phi_coeff';
        
        %chain rule for derivatives
        derivative_1 = N(indexLocalEnrichedNodes(nodalIndex)) * Phi_iModeContinuousDer;
        derivative_2 =  B(indexLocalEnrichedNodes(nodalIndex)) * Phi_iMode;
        
        %PhiDerivative modal basis derivative
        PhiDerivative = [PhiDerivative, derivative_1 + derivative_2];
    end
end
end

