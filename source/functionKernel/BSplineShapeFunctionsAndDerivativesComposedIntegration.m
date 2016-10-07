function [ N, B ] = BSplineShapeFunctionsAndDerivativesComposedIntegration( x,...
    integrationSubDomainIndex, coefficients, leftKnot, rightKnot,...
    leftIntegrationDomain, rightIntegrationDomain, problem )
%BSPLINESHAPEFUNCTIONSANDDERIVATIVESCOMPOSEDINTEGRATION evaluate BSpline
%functions and their derivatives at x. 
%Input:
% x = integration point local coordinate
%integrationSubDomainIndex = index of the integration domain
%coefficients = coefficients
%leftKnot = left knot of the knot span
%rightKnot = right knot of the knot span
%leftIntegrationDomain = left node of the integration subDomain
%rightIntegrationDomain = right node of the integration subDomain
%problem = base mesh problem struct
%Output:
%N = Bspline shape functions
%B = BSplines shape function derivatives 

if problem.p == 1
    %% Matrix of shape functions
    
    localCoordinates = mapGlobalToLocal(x, leftIntegrationDomain, rightIntegrationDomain);

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
    
else
    
    m = length(coefficients);
    evaluationBSplinesPoint = linspace(leftKnot, rightKnot, m);
    localCoordinates = mapGlobalToLocal(x, leftIntegrationDomain, rightIntegrationDomain);
    
    N = [];
    B = [];
    
    [NSpline1, ~] = BsplinesShapeFunctionsAndDerivatives(evaluationBSplinesPoint(integrationSubDomainIndex),...
        problem.p, problem.knotVector );
    [NSpline2, ~] = BsplinesShapeFunctionsAndDerivatives(evaluationBSplinesPoint(integrationSubDomainIndex + 1),...
        problem.p, problem.knotVector );

    N_coeff1 = NSpline1(end-problem.p:end);
    N_coeff2 = NSpline2(end-problem.p:end);

    %% Matrix of shape functions
    
    N_subElement =  0.5 * [(1 - localCoordinates), (1 + localCoordinates)];
    for i =1:problem.p+1
        N_vector = N_subElement(1) * N_coeff1(i) + N_subElement(2) * N_coeff2(i);
        N = [N, N_vector];
    end
    
    %% Matrix of shape functions derivatives
    
    B_subElement =  0.5 * [-1, 1];  

    for i =1:problem.p+1
        B_vector = B_subElement(1) * N_coeff1(i) + B_subElement(2) * N_coeff2(i);
        B = [B, B_vector];
    end
end

end