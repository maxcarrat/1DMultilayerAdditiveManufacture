function [ numericalSolutions ] = evaluateNumericalResultsIGA( x, t, problem,...
    coefficients, layer, numberOfLayers, derivative )
% numericalSolutions = evaluateNumerialResultsIGA(postProcessingCoords, problem, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% problem = struct that defines the boundary value problem
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% derivative = index of the deriative of the element numerical derivative to be evaluated

numericalSolutions=zeros(size(x));
Xi1 = problem.knotVector( 1 + problem.p );
Xi2 = problem.knotVector( 2 + problem.p );

x = linspace(0, 1, length(x) / numberOfLayers * layer );
numericalSolutions(x>=Xi1 & x<=Xi2) = element_num_sol( x(x>=Xi1 & x<=Xi2), t, problem, 1, coefficients, derivative);

for e=2:layer-1
    Xi1 = problem.knotVector( e + problem.p );
    Xi2 = problem.knotVector( e + 1 + problem.p );
      
numericalSolutions(x>Xi1 & x<=Xi2) = element_num_sol( x(x>Xi1 & x<=Xi2), t, problem, e, coefficients, derivative);

end

for e=layer+(problem.p-1):size(problem.LM, 1)
    Xi1 = problem.knotVector( e + problem.p );
    Xi2 = problem.knotVector( e + 1 + problem.p );
      
numericalSolutions(x>Xi1 & x<=Xi2) = element_num_sol( x(x>Xi1 & x<=Xi2), t, problem, e, coefficients, derivative);

end

end

function r = element_num_sol(x, t, problem, element, coefficients, derivative)
% r = ELEMENT_NUM_SOL(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated

numberOfProjectionPoints = length(x);
m = length(problem.knotVector);

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(x), m - 1 - problem.p);
B = zeros(length(x), m - 1 - problem.p);
JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

for k=1:length(x)
end

if derivative == 0
    
    for i=1:length(x)
        [N(k,:), ~] = BsplinesShapeFunctionsAndDerivatives(x(k), problem.p, problem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
    end
    
    r = projectionOperator(:,problem.LM(element,:)) * coefficients(problem.LM(element,:)) ;

else
    
    for i=1:length(x)
        [N(k,:), B(k,:)] = BsplinesShapeFunctionsAndDerivatives(x(k), problem.p, problem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
        projectionOperator_der(i,1:size(B,2)) = B(i,:);
        JacobianX_Xi(i) = B(i,problem.LM(element, :)) *...
            problem.coords(problem.LM(element, :))';
        inverseJacobianX_Xi(i) = 1 / JacobianX_Xi(i);
    end
   
    r = (projectionOperator_der(:,problem.LM(element,:)) * coefficients(problem.LM(element,:))) .*...
        inverseJacobianX_Xi .*  ...
        problem.k(mapParametricToGlobal(x, problem), t, projectionOperator(:,problem.LM(element,:)) * coefficients(problem.LM(element,:))) ;

end

end