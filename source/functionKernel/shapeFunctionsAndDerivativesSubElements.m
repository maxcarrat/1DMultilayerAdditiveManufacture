function [ N, B ] = shapeFunctionsAndDerivativesSubElements( x,...
    integrationSubDomainIndex, coefficients, integrationDomain )
%SHAPEFUNCTIONSANDDERIVATIVESSUBELEMENTS  Evaluate the ahpe functions and their derivatives 
% on the integration sub-elements 
%   x = integration point local coordinate
%   coefficients = nodal values of shape function at the element level
%   integrationSubDomainIndex = index of the integration sub-domain

%% Matrix of shape functions
leftElementNode = integrationDomain(integrationSubDomainIndex);
rightElementNode = integrationDomain(integrationSubDomainIndex + 1);
localCoordinates = mapGlobalToLocal(x, leftElementNode, rightElementNode);

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

subDomainBackwardMapping = 2 / (rightElementNode - leftElementNode);

B = [B_subElement * ( 1 - B_coeff)', B_subElement *  B_coeff'] * subDomainBackwardMapping;


end

