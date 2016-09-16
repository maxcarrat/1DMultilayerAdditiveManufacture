function [ N, B ] = shapeFunctionsAndDerivativesIGASubElements( localCoordinates,...
    integrationSubDomainIndex, coefficients, leftKnot, rightKnot, problem )
%SHAPEFUNCTIONSANDDERIVATIVESIGASUBELEMENTS  Evaluate the shape functions and their derivatives 
% on the integration sub-elements 
%   localCoordinate = integration point local coordinate
%   coefficients = nodal values of shape function at the element level
%   integrationSubDomainIndex = index of the integration sub-domain
%   leftKnot = knot on the elemnt left end side
%   rightKnot = knot on the elemnt right end side


if problem.p == 1
    %% Matrix of shape functions

    N_subElement =  0.5 * [(1 - localCoordinates), (1 + localCoordinates)];
    N1 = coefficients(integrationSubDomainIndex);
    N2 = coefficients(integrationSubDomainIndex+1);
    N_coeff = [N1, N2];
    
    N = [N_subElement * ( 1 - N_coeff)', N_subElement * N_coeff'];
    
    
    %% Matrix of shape functions derivatives
    
    inverseJacobianDeterminantParentToParametric = 2/(rightKnot - leftKnot);
    
    B_subElement =  0.5 * [-1, 1];
    B1 = coefficients(integrationSubDomainIndex);
    B2 = coefficients(integrationSubDomainIndex+1);
    B_coeff = [B1, B2];
    
    B = [B_subElement * ( 1 - B_coeff)', B_subElement *  B_coeff'] * inverseJacobianDeterminantParentToParametric;
    
else
    
    m = length(coefficients);
    evaluationBSplinesPoint = linspace(leftKnot, rightKnot, m);
    
    N = [];
    B = [];
    
    [NSpline1, BSpline1] = BsplinesShapeFunctionsAndDerivatives(evaluationBSplinesPoint(integrationSubDomainIndex),...
        problem.p, problem.knotVector );
    [NSpline2, BSpline2] = BsplinesShapeFunctionsAndDerivatives(evaluationBSplinesPoint(integrationSubDomainIndex + 1),...
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
    
    inverseJacobianDeterminantParentToParametric = 2/(rightKnot - leftKnot);
    
    B_subElement =  0.5 * [-1, 1]; 
    B_coeff1 = NSpline1(end-problem.p:end);
    B_coeff2 = NSpline2(end-problem.p:end);    

    for i =1:problem.p+1
        B_vector = B_subElement(1) * B_coeff1(i) + B_subElement(2) * B_coeff2(i);
        B = [B, B_vector * inverseJacobianDeterminantParentToParametric ];
    end
end

end


