function [ N, B ] = shapeFunctionsAndDerivativesSubElements( localCoordinates,...
    integrationSubDomainIndex, coefficients )
%SHAPEFUNCTIONSANDDERIVATIVESSUBELEMENTS  Evaluate the ahpe functions and their derivatives 
% on the integration sub-elements 
%   localCoordinate = integration point local coordinate
%   coefficients = nodal values of shape function at the element level
%   integrationSubDomainIndex = index of the integration sub-domain

%% Matrix of shape functions

N_subElement =  0.5 * [(1 - localCoordinates), (1 + localCoordinates)];
N1 = coefficients(integrationSubDomainIndex);
N2 = coefficients(integrationSubDomainIndex+1);
N_coeff = [N1, N2];

N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];


%% Matrix of shape functions derivatives

B_subElement =  0.5 * [-1, 1];
B1 = coefficients(integrationSubDomainIndex);
B2 = coefficients(integrationSubDomainIndex+1);
B_coeff = [B1, B2];

B = [B_subElement * ( 1 - B_coeff)', B_subElement *  B_coeff'];


end

