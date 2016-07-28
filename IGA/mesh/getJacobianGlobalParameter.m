function [ J ] = getJacobianGlobalParameter( x, problem , derivativesOfBSpline, e)
%GETJACOBIANGLOBALPARAMETER

i = findspan(length(problem.coords)-1, problem.p, x, problem.knotVector);
J = 0.0;

for j=1:problem.p+1
    
    J = J + (derivativesOfBSpline((e-1)+j) * problem.coords(i - problem.p + j))^2;

end

J = sqrt(J);

end

