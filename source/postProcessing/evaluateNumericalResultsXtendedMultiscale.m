function [ numericalSolutions ] = evaluateNumericalResultsXtendedMultiscale(  x, t, baseProblem,...
    initialOverlayProblem, overlayProblem, baseCoefficients, overlayCoefficients,...
    layer, numberOfLayers, derivative  )
%EVALUATENUMERICALRESULTSXTENDEDMULTISCALE   evaluates the numerical solution of
%the problem.
%Input:
%x = coordsinates to post process
%t = actual time
%baseProblem = struct that defines the boundary value problem on the base
%mesh
%initialOverlayProblem = struct of the overlay hFEM problem
%overlayProblem = struct that defines the boundary value problem on the
%overlay mesh using X-PODFEM
%baseCoefficients = coefficients of the basis function on the base mesh
%overlayCoefficients = coefficients of the basis function on the overlay
%mesh
%layer = layer of the AM process
%numberOfLayers = total number of layers
%derivative = index of the deriative of the element numerical derivative to
%be evaluated
%Output:
%numericalSolutions = solution value at the post-processing points

%% Initialize variables
numericalSolutions=zeros(size(x));
Xi1 = baseProblem.knotVector( 1 + baseProblem.p );
Xi2 = baseProblem.knotVector( 2 + baseProblem.p );

xParametric = linspace(0, 1, length(x) / numberOfLayers * layer );
numericalSolutions(xParametric>=Xi1 & xParametric<=Xi2) = numericalSolutionBaseMesh( xParametric(xParametric>=Xi1 & xParametric<=Xi2), t, baseProblem,...
    baseCoefficients, 1, derivative);

%loop over non-overlapped elements
for e=2:layer-1
    Xi1 = baseProblem.knotVector( e + baseProblem.p );
    Xi2 = baseProblem.knotVector( e + 1 + baseProblem.p );
    
    numericalSolutions(xParametric>Xi1 & xParametric<=Xi2) = numericalSolutionBaseMesh( xParametric(xParametric>Xi1 & xParametric<=Xi2), t, baseProblem,...
        baseCoefficients, e, derivative);
end

%loop over overlapped elements
X1 = overlayProblem.coords(1);
X2 = overlayProblem.coords(2);
numericalSolutions(x>=X1 & x<=X2) = numericalSolutionOverlayMesh( x(x>=X1 & x<=X2), t, baseProblem,...
    overlayProblem, baseCoefficients, initialOverlayProblem, overlayCoefficients, 1, derivative);

for e=2:overlayProblem.N
    X1 = overlayProblem.coords(e);
    X2 =overlayProblem.coords(e+1);
    
    numericalSolutions(x>X1 & x<=X2) = numericalSolutionOverlayMesh( x(x>X1 & x<=X2), t, baseProblem,...
        overlayProblem, baseCoefficients, initialOverlayProblem, overlayCoefficients, e, derivative);
end

end

function r = numericalSolutionBaseMesh(x, t, baseProblem, ...
    baseCoefficients, element, derivative)
%NUMERICALSOLUTIONBASEMESH evaluates the numerical solution associated with
%a single specific element of the base mesh
%Input:
%x = points where the element numerical solution has to be evaluated
%t = actual time
%baseProblem = struct that defines the boundary value problem on the base
%mesh
%baseCoefficients = coefficients of the basis function on the base mesh
%element = index of the element where to evaluate the element numerical
%solution
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

%% Initialize variables

numberOfProjectionPoints = length(x);
m = length(baseProblem.knotVector);

projectionOperator = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));

N = zeros(length(x), m - 1 - baseProblem.p);
B = zeros(length(x), m - 1 - baseProblem.p);
JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

%% Interpolate solution coefficients by means of BSplines and derivatives

%evaluate temperature
if derivative == 0
    for i=1:length(x)
        [N(i,:), ~] = BsplinesShapeFunctionsAndDerivatives(x(i), baseProblem.p, baseProblem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
    end
    
    r = projectionOperator(:,baseProblem.LM(element,:)) * baseCoefficients(baseProblem.LM(element,:)) ;
    
    %evaluate heat fluxes
else
    for i=1:length(x)
        [N(i,:), B(i,:)] = BsplinesShapeFunctionsAndDerivatives(x(i), baseProblem.p, baseProblem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
        projectionOperator_der(i,1:size(B,2)) = B(i,:);
        
        JacobianX_Xi(i) = baseProblem.coords(end) - baseProblem.coords(1);
        inverseJacobianX_Xi(i) = 1 / JacobianX_Xi(i);
    end
    
    r = (projectionOperator_der(:,baseProblem.LM(element,:)) * baseCoefficients(baseProblem.LM(element,:))) .*...
        inverseJacobianX_Xi .*  ...
        baseProblem.k(x, t, projectionOperator(:,baseProblem.LM(element,:)) *...
        baseCoefficients(baseProblem.LM(element,:)));
end

end

function r = numericalSolutionOverlayMesh(x, t, baseProblem, overlayProblem, ...
    baseCoefficients, initialOverlayProblem, overlayCoefficients, element, derivative)
%NUMERICALSOLUTIONOVERLAYMESH evaluates the numerical solution associated with
%a single specific element of the overlay mesh
%Input:
%x = points where the element numerical solution has to be evaluated
%t = actual time
%baseProblem = struct that defines the boundary value problem on the base
%mesh
%overlayProblem = struct that defines the boundary value problem on the
%overlay mesh
%baseCoefficients = coefficients of the basis function on the base mesh
%initialOverlayProblem = struct of the overlay hFEM problem
%overlayCoefficients = coefficients of the basis function on the overlay
%mesh
%element = index of the element where to evaluate the element numerical
%solution
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

%% Initialize variables ---------------------------------------------------
%Base Mesh
numberOfProjectionPoints = length(x);
m = length(baseProblem.knotVector);

projectionOperatorSplines = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
projectionOperatorSplines_der = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));

JacobianX_Xi = zeros(length(x), 1);

numberOfProjectionPoints = length(x);

%Overlay Mesh
projectionOperator_der = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2));

parametricCoords = mapGlobalToParametric(x, baseProblem.coords(1), baseProblem.coords(end));

%evaluate temperature

%get the local evaluation point
localCoords = mapGlobalToLocal(x, overlayProblem.coords(element), overlayProblem.coords(element+1));
refinedDofs = initialOverlayProblem.N+1;
integrationDomain = linspace(-1, +1, ceil(refinedDofs/overlayProblem.XN));

%overlay mesh element boundaries
X1 = overlayProblem.coords(element);
X2 = overlayProblem.coords(element +1);
elementGlobalCoords = [X1, X2];

inverseJacobianLocalToGlobal = overlayProblem.B_map(X1,X2);
    
elementEnrichedIndex = element - (overlayProblem.N - overlayProblem.XN);

if elementEnrichedIndex == 1
    indexLocalEnrichedNodes = 2;
else
    indexLocalEnrichedNodes = [1, 2];
end

modalDofs = length(indexLocalEnrichedNodes)*overlayProblem.modes;

projectionOperator = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2) + modalDofs);

r = zeros(1, numberOfProjectionPoints);
r_base = zeros(1, numberOfProjectionPoints);
r_over = zeros(1, numberOfProjectionPoints);

enrichedSubElements = refinedDofs/overlayProblem.XN;

PODCoefficients = overlayProblem.reductionOperator(...
    (elementEnrichedIndex-1)*floor(enrichedSubElements)+1:(elementEnrichedIndex-1)*floor(enrichedSubElements)...
    + ceil(enrichedSubElements),:);

%for each post-processing point...
for i=1:numberOfProjectionPoints
    
    %...loop over sub-domains
    for integrationSubDomainIndex = 1 : ceil(enrichedSubElements) - 1
        
        %if the GP is inside the sub-domain
        if localCoords(i) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                localCoords(i) > integrationDomain(integrationSubDomainIndex)
            
            if derivative == 0
                
                [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(parametricCoords(i),...
                    baseProblem.p, baseProblem.knotVector);
                projectionOperatorSplines(i,1:length(NIga)) = NIga(:);
                [N, ~] =  shapeFunctionsAndDerivatives(localCoords(i));
                [F, ~] = PODModesAndDerivativesMultiscale(localCoords(i), elementGlobalCoords, overlayProblem.modes,...
                    PODCoefficients, integrationDomain, integrationSubDomainIndex, indexLocalEnrichedNodes);
                projectionOperator(i,1:size(N,2)+size(F,2)) = [N, F];
                
                
                r_base(i) = projectionOperatorSplines(i,:) * baseCoefficients(:) ;
                r_over(i) = projectionOperator(i,:) * overlayCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2));

                %the final solution is the sum of the base and overlay mesh solutions
                r(i) = r_base(i) + r_over(i);

                %evaluate heat fluxes
            else
                
                JacobianX_Xi(i) = baseProblem.coords(end) - baseProblem.coords(1);
                inverseJacobianX_Xi = 1 / JacobianX_Xi(i);
                
                [NIga, BIga] =  BsplinesShapeFunctionsAndDerivatives(parametricCoords(i),...
                    baseProblem.p, baseProblem.knotVector);
                
                projectionOperatorSplines(i,1:length(NIga)) = NIga(:);
                projectionOperatorSplines_der(i,1:length(BIga)) = BIga(:) .* inverseJacobianX_Xi;
                
                
                [N, B] =  shapeFunctionsAndDerivatives(localCoords(i));
                B = B * inverseJacobianLocalToGlobal;
                
                [F, G] = PODModesAndDerivativesMultiscale(localCoords(i), elementGlobalCoords, overlayProblem.modes,...
                    PODCoefficients, integrationDomain, integrationSubDomainIndex, indexLocalEnrichedNodes);
                
                projectionOperator(i,1:size(N,2)+size(F,2)) = [N, F];
                projectionOperator_der(i,1:size(B,2)+size(G,2)) = [B, G];
                
                temperature = projectionOperator(i,:) * overlayCoefficients(...
                    overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)) + projectionOperatorSplines(i,:) *...
                    baseCoefficients(:);
                
                r_base(i) = (projectionOperatorSplines_der(i,:) * baseCoefficients(:)) .* ...
                    baseProblem.k(x(i), t, temperature);
                
                r_over(i) = (projectionOperator_der(i,:) * overlayCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2))) .*...
                    baseProblem.k(x(i), t, temperature);
               
                %the final solution is the sum of the base and overlay mesh solutions
                r(i) = r_base(i) + r_over(i);
            end
            break;
        end
        
        
    end %end loop over subDomains
    
    
    
end %end loop over post-processing points


end


