function [ numericalSolutions ] = evaluateNumericalResultsEnriched( x, activeElementCoords, problem, Xcoefficients, coefficients, derivative )
% numericalSolutions = evaluateNumerialResults(postProcessingCoords, problem, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% problem = struct that defines the boundary value problem
% Xcoefficients = coefficients of the enriched solution
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% derivative = index of the deriative of the element numerical derivative to be evaluated

coords = problem.coords;
numericalSolutions=zeros(size(x));
N = size(problem.LM, 1);
X1 = coords(1);
X2 = coords(2);

if(X1>=activeElementCoords(end-1)) && (X2 <= activeElementCoords(end))
    numericalSolutions(x>=X1 & x<=X2) = element_num_sol_enriched(x(x>=X1 & x<=X2), coords, problem, 1, Xcoefficients, derivative);
else
    numericalSolutions(x>=X1 & x<=X2) = element_num_sol_linear(x(x>=X1 & x<=X2), coords, problem, 1, coefficients, derivative);
end

for e=2:N
    X1 = coords(e);
    X2 = coords(e+1);
    
    if(e == N)
        numericalSolutions(x>X1 & x<=X2) = element_num_sol_enriched(x(x>X1 & x<=X2), coords, problem, e, Xcoefficients, derivative);
    else
        numericalSolutions(x>X1 & x<=X2) = element_num_sol_linear(x(x>X1 & x<=X2), coords, problem, e, coefficients, derivative);
    end
end

end

function r = element_num_sol_linear(x, coords, problem, element, coefficients, derivative)
% r = ELEMENT_NUM_SOL_LINEAR(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated

X1 = coords(element);
X2 = coords(element+1);
ldof = 2;

r = zeros(size(x));
for i=1:ldof
    r=r+coefficients(problem.LM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, derivative);
end

r = r .* (2/(X2-X1)) ^ derivative;

end

function r = element_num_sol_enriched(x, coords, problem, element, coefficients, derivative)
% r = ELEMENT_NUM_SOL_ENRICHED(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated

X1 = coords(element);
X2 = coords(element+1);
ldof = 2;

enrichedElementCoords = linspace(X1, X2, 2^(problem.refinementDepth)+1);
r = zeros(size(x));

if derivative==0
    
    for i=1:ldof
        r=r+coefficients(problem.LM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, 0);
    end
    for iMode=1:problem.modes
        for i =1:ldof
            r=r+coefficients(problem.N+1+(i-1)*problem.modes+iMode).*...
                problem.xFEMBasis_fun(mapGlobalToLocal(x, X1, X2), i, iMode, 0, 0,...
                problem, enrichedElementCoords);
        end
    end
    
else
    
    for i=1:ldof
        r=r+coefficients(problem.LM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, 1);
    end
    for iMode=1:problem.modes
        for i =1:ldof
            r=r+coefficients(problem.N+1+(i-1)*problem.modes+iMode).*...
                (problem.xFEMBasis_fun(mapGlobalToLocal(x, X1, X2), i, iMode, 1, 0, problem, enrichedElementCoords) + ...
                problem.xFEMBasis_fun(mapGlobalToLocal(x, X1, X2), i, iMode, 0, 1, problem, enrichedElementCoords));
        end
    end
    
end

% r = r .* (2/(X2-X1)) ^ derivative;

end



