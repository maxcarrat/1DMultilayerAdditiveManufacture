function [ numericalSolutions ] = evaluateNumericalResultsMultiscale( x, t, baseProblem,...
    overlayProblem, baseCoefficients, overlayCoefficients,...
    layer, numberOfLayers, derivative )
%EVALUATENUMERICALRESULTSMULTISCALE  evaluates the numerical solution of
%the problem.
%Input:
%x = coordsinates to post process
%t = actual time
%baseProblem = struct that defines the boundary value problem on the base
%mesh
%overlayProblem = struct that defines the boundary value problem on the
%overlay mesh
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
X1 = overlayProblem.coords( 1 );
X2 = overlayProblem.coords( 2 );

numericalSolutions(x>=X1 & x<=X2) = numericalSolutionOverlayMesh( x(x>=X1 & x<=X2), t, baseProblem,...
    overlayProblem, baseCoefficients, overlayCoefficients, 1, derivative);

for e=2:size(overlayProblem.LM, 1)
    X1 = overlayProblem.coords( e );
    X2 = overlayProblem.coords( e + 1 );

    numericalSolutions(x>X1 & x<=X2) = numericalSolutionOverlayMesh( x(x>X1 & x<=X2), t, baseProblem,...
        overlayProblem, baseCoefficients, overlayCoefficients, e, derivative);
end

end

function r = numericalSolutionBaseMesh(x, t, baseProblem, ...
    baseCoefficients, element, derivative)
%NUMERICALSOLUTIONBASEMESH evaluates the numerical solution associated with 
%a single specific element of the base mesh
%Input:
%x = parametric coords of the points where the element numerical solution has to be evaluated
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

JacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
inverseJacobianX_Xi = 1 / JacobianX_Xi;

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
        
        
    end
    
    temp = projectionOperator(:,baseProblem.LM(element,:)) * baseCoefficients(baseProblem.LM(element,:));
    
    
    r = (projectionOperator_der(:,baseProblem.LM(element,:)) * baseCoefficients(baseProblem.LM(element,:))) .*...
        inverseJacobianX_Xi .* baseProblem.k(mapParametricToGlobal(x, baseProblem), t, temp);
end

end

function r = numericalSolutionOverlayMesh(x, t, baseProblem, overlayProblem, ...
    baseCoefficients, overlayCoefficients, element, derivative)
%NUMERICALSOLUTIONOVERLAYMESH evaluates the numerical solution associated with 
%a single specific element of the overlay mesh
%Input:
%x = global coords of the points where the element numerical solution has to be evaluated
%t = actual time
%baseProblem = struct that defines the boundary value problem on the base
%mesh
%overlayProblem = struct that defines the boundary value problem on the
%overlay mesh
%baseCoefficients = coefficients of the basis function on the base mesh
%overlayCoefficients = coefficients of the basis function on the overlay
%mesh%element = index of the element where to evaluate the element numerical
%solution
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

%% Initialize variables ---------------------------------------------------
%Base Mesh
numberOfProjectionPoints = length(x);

projectionOperatorSplines = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
projectionOperatorSplines_der = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));

Nspline = zeros(length(x), size(baseProblem.LM, 2));
Bspline = zeros(length(x), size(baseProblem.LM, 2));


JacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
inverseJacobianX_Xi = 1 / JacobianX_Xi;

%Overlay Mesh
X1 = overlayProblem.coords(element);
X2 = overlayProblem.coords(element+1);

numberOfProjectionPoints = length(x);

projectionOperator = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2));

N = zeros(length(x), size(overlayProblem.LM, 2));
B = zeros(length(x), size(overlayProblem.LM, 2));

globalcoords = x;
localcoords = mapGlobalToLocal( globalcoords, X1, X2);

%evaluate temperature
if derivative == 0
    for i=1:length(x)
        [N, ~] =  shapeFunctionsAndDerivatives(localcoords(i));
        [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(globalcoords(i),...
            baseProblem.coords(1),baseProblem.coords(end)),...
            baseProblem.p, baseProblem.knotVector);
        projectionOperatorSplines(i,1:size(Nspline,2)) = NIga(baseProblem.LM(end,:));
        projectionOperator(i,1:size(N,2)) = N;
    end
    
    r_base = projectionOperatorSplines(:,:) * baseCoefficients(baseProblem.LM(end,:)) ;
    r_over = projectionOperator * overlayCoefficients(overlayProblem.LM(element,:)) ;
    
    %evaluate heat fluxes
else
    for i=1:length(x)
        [N, B] = shapeFunctionsAndDerivatives(localcoords(i));
        [NIga, BIga] =  BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(globalcoords(i),...
            baseProblem.coords(1),baseProblem.coords(end)),...
            baseProblem.p, baseProblem.knotVector);
        projectionOperatorSplines(i,1:size(Nspline,2)) = NIga(baseProblem.LM(end,:));
        projectionOperatorSplines_der(i,1:size(Bspline,2)) = BIga(baseProblem.LM(end,:));
        projectionOperator(i,1:size(N,2)) = N;
        projectionOperator_der(i,1:size(B,2)) = B;
    end
    
    temp = projectionOperatorSplines(:,:) *...
        baseCoefficients(baseProblem.LM(end,:)) +  projectionOperator * overlayCoefficients(overlayProblem.LM(element,:));
    
    r_base = (projectionOperatorSplines_der(:,:) * baseCoefficients(baseProblem.LM(end,:))) .*...
        inverseJacobianX_Xi^derivative .* baseProblem.k(globalcoords, t, temp);
    
    r_over = (projectionOperator_der * overlayCoefficients(overlayProblem.LM(element,:))).*...
        overlayProblem.B_map(X1, X2)^ derivative .* overlayProblem.k(globalcoords, t, temp); 
end

%the final solution is the sum of the base and overlay mesh solutions
r = r_base + r_over;

end

