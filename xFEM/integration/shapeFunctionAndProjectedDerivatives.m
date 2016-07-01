function [ N, B ] = shapeFunctionAndProjectedDerivatives( localCoordinates,...
                        coefficients, shapeFunctCoeff, integrationSubDomainIndex )
%SHAPEFUNCTIONANDPROJECTEDDERIVATIVES Summary of this function goes here
%   Detailed explanation goes here

N = 0.5 * [(1-localCoordinates), (1+localCoordinates)];

%% Matrix of shape functions derivatives

Xi1 = coefficients(integrationSubDomainIndex);
Xi2 = coefficients(integrationSubDomainIndex + 1);
localSubElementCoords = mapGlobalToLocal( localCoordinates, Xi1, Xi2);

N_subElement =  0.5 * [(1 - localSubElementCoords)', (1 + localSubElementCoords)'];
B_subElement =  0.5 * [ -1, 1];

%POD coefficients vector
q_1 = shapeFunctCoeff(integrationSubDomainIndex);
q_2 = shapeFunctCoeff(integrationSubDomainIndex+1);
q_coeff = [q_1, q_2];

der = B_subElement * ( q_coeff )';

% Galerkin projection
[rGP, wGP] = gaussPoints( 2 );
M = zeros(2,2);
b = zeros(2,1);

for iGP = 1:numel(rGP)
    [n, ~] = shapeFunctionsAndDerivatives(rGP(iGP));
    M = n' * n * wGP(iGP) + M;
    b = n' * der * wGP(iGP) + b;
end

continuousDer = M\b;

B = [N_subElement(1) * continuousDer(1), N_subElement(2) * continuousDer(2)];

end