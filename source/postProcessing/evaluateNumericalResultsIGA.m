function [ numericalSolutions ] = evaluateNumericalResultsIGA( x, t, problem,...
    coefficients, layer, numberOfLayers, derivative )
%EVALUATENUMERICALRESULTSIGA evaluates the numerical solution
%Input:
%x = coordinates to post process
%t = actual time
%problem = struct that defines the boundary value problem
%coefficients = coefficients of the basis function obtained by solving the
%mass matrix-load vector system of equations
%layer = layer of the AM process
%numberOfLayers = total number of layers
%derivative = index of the deriative of the element numerical derivative to
%be evaluated
%Output:
%numericalSolutions = solution value at the post-processing points

%% Initialize variables
numericalSolutions=zeros(size(x));
numberOfElements = problem.N;

Xi1 = problem.knotVector( 1 + problem.p );
Xi2 = problem.knotVector( 2 + problem.p );

X1 = mapParametricToGlobal(Xi1, problem);
X2 = mapParametricToGlobal(Xi2, problem);

numericalSolutions(x>=X1 & x<=X2) = element_num_sol( x(x>=X1 & x<=X2), t, problem, 1,...
    coefficients, derivative, X1, X2 );

%loop over non-enriched elements
for e=2:numberOfElements
    Xi1 = problem.knotVector( e + problem.p );
    Xi2 = problem.knotVector( e + 1 + problem.p );
    
    X1 = mapParametricToGlobal(Xi1, problem);
    X2 = mapParametricToGlobal(Xi2, problem);
    
    numericalSolutions(x>X1 & x<=X2) = element_num_sol( x(x>X1 & x<=X2), t, problem, e,...
        coefficients, derivative, X1, X2 );
end

% %loop over enriched elements
% for e=layer+(problem.p-1):size(problem.LM, 1)
%     Xi1 = problem.knotVector( e + problem.p );
%     Xi2 = problem.knotVector( e + 1 + problem.p );
%       
% numericalSolutions(x>Xi1 & x<=Xi2) = element_num_sol( x(x>Xi1 & x<=Xi2), t, problem, e, coefficients, derivative);
% end

end

function r = element_num_sol(x, t, problem, element, coefficients, derivative, X1, X2 )
%ELEMENT_NUM_SOL evaluates the numerical solution associated with a single specific element
%Input:
% x = points where the element numerical solution has to be evaluated
%t = actual time
%problem = struct that defines the boundary value problem
%element = index of the element where to evaluate the element numerical
%solution
% coefficients = coefficients of the basis function obtained by solving the
%mass matrix-load vector system of equations
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

%% Initialize variables
numberOfProjectionPoints = length(x);
m = length(problem.knotVector);

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(x), m - 1 - problem.p);
B = zeros(length(x), m - 1 - problem.p);

JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

parametricCoords = mapGlobalToParametric(x, problem.coords(1), problem.coords(end));

%% Interpolate solution coefficients by means of BSplines and derivatives

%evaluate temperature
if derivative == 0
    for i=1:length(x)
        [N(i,:), ~] = BsplinesShapeFunctionsAndDerivatives(parametricCoords(i), problem.p, problem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
    end
    
    r = projectionOperator(:,problem.LM(element,:)) * coefficients(problem.LM(element,:)) ;
   
%evaluate heat fluxes
else
    for i=1:length(x)
        [N(i,:), B(i,:)] = BsplinesShapeFunctionsAndDerivatives(parametricCoords(i), problem.p, problem.knotVector);
        projectionOperator(i,1:size(N,2)) = N(i,:);
        projectionOperator_der(i,1:size(B,2)) = B(i,:);
        
        JacobianX_Xi(i) = B(i,problem.LM(element, :)) *...
            problem.coords(problem.LM(element, :))';
        
        %         detJacobianX_Xi = norm(JacobianX_Xi);
        detJacobianX_Xi = problem.coords(end) - problem.coords(1);
        
        inverseJacobianX_Xi(i) = 1 / detJacobianX_Xi;
    end
    
    r = (projectionOperator_der(:,problem.LM(element,:)) * coefficients(problem.LM(element,:))) .*...
        inverseJacobianX_Xi .*  ...
        problem.k(mapLocalToGlobal(x, X1, X2), t, projectionOperator(:,problem.LM(element,:)) *...
        coefficients(problem.LM(element,:)));
end

end