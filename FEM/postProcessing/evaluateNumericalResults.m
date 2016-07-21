function [ numericalSolutions ] = evaluateNumericalResults( x, t, problem, coefficients, derivative )
% numericalSolutions = evaluateNumerialResults(postProcessingCoords, problem, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% problem = struct that defines the boundary value problem
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% derivative = index of the deriative of the element numerical derivative to be evaluated    

    coords = problem.coords;
    numericalSolutions=zeros(size(x));
    X1 = coords(1);
    X2 = coords(2);

    numericalSolutions(x>=X1 & x<=X2) = element_num_sol(x(x>=X1 & x<=X2), coords, t, problem, 1, coefficients, derivative);
    
    for e=2:size(problem.LM, 1)
        X1 = coords(e);
        X2 = coords(e+1);
        
        numericalSolutions(x>X1 & x<=X2) = element_num_sol(x(x>X1 & x<=X2), coords, t, problem, e, coefficients, derivative);
    end
end

function r = element_num_sol(x, coords, t, problem, element, coefficients, derivative)
% r = ELEMENT_NUM_SOL(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated    

X1 = coords(element);
X2 = coords(element+1);

numberOfProjectionPoints = length(x);

 projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
 projectionOperator_der = zeros(numberOfProjectionPoints, size(problem.LM, 2));
    
    N = zeros(length(x), 2);
    B = zeros(length(x), 2);
    
    localCoords = mapGlobalToLocal( x, X1, X2);
    
    for k=1:length(x)
        [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    end
    
    if derivative == 0
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
        r = projectionOperator * coefficients(problem.LM(element,:)) ;
    else
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
            projectionOperator_der(i,1:size(B,2)) = B(i,:);
        end

        r = (projectionOperator_der * coefficients(problem.LM(element,:))).* (2/(X2-X1)) ^ derivative .*...
                problem.k(x, t, projectionOperator * coefficients(problem.LM(element,:))) ;
    end
    
%     r = projectionOperator * coefficients(problem.LM(element,:)) ;
    
    %     for i=1:ldof
    %         r=r+coefficients(problem.LM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, derivative);
    %     end
    %

    
%     r = r .* (2/(X2-X1)) ^ derivative + derivative * problem.k(mapGlobalToLocal(x, X1, X2), t, projectionOperator * coefficients(problem.LM(element,:))) ;

end