function [ numericalSolutions ] = evaluateXFEMNumericalResults( x, problem, coefficients, derivative )
% numericalSolutions = evaluateXFEMNumerialResults(postProcessingCoords, problem, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% problem = struct that defines the boundary value problem
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% derivative = index of the deriative of the element numerical derivative to be evaluated    

    coords = problem.coords;
    numericalSolutions=zeros(size(x));
    X1 = coords(1);
    X2 = coords(2);

    numericalSolutions(x>=X1 & x<=X2) = element_num_sol(x(x>=X1 & x<=X2), coords, problem, 1, coefficients, derivative);
    
    for e=2:problem.N
        X1 = coords(e);
        X2 = coords(e+1);
        
        numericalSolutions(x>X1 & x<=X2) = element_num_sol(x(x>X1 & x<=X2), coords, problem, e, coefficients, derivative);
    end
end

function r = element_num_sol(x, coords, problem, element, coefficients, derivative)
% r = ELEMENT_NUM_SOL(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated

X1 = coords(element);
X2 = coords(element+1);

r = zeros(size(x));

for i=1:2
    
    r=r+coefficients(problem.rbLM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, derivative);
end

for iMode = 3:2*problem.modes+2
    
    r=r+coefficients(problem.rbLM(element,i)).* problem.xFEMBasis_fun(mapGlobalToLocal(x, X1, X2), i, iMode, derivative,...
        derivative, problem, problem.coords);
end


r = r .* (2/(X2-X1)) ^ derivative;

end

