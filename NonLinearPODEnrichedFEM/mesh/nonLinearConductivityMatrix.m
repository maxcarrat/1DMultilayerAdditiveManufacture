function [ KE ] = nonLinearConductivityMatrix( element, problem, actualTemperature )
%NONLINEARCONDUCTIVITYMATRIX computes the conductivity matrix of the
%element as a function of the actual temperature
%   element = actual element of the mesh
%   problem = non-linear transient heat problem struct
%   actualTemperature = temperature valueas at the current time step

X1 = problem.coords(element);
X2 = problem.coords(element+1);

ldof = 2;
KE = zeros(ldof, ldof);

for i=1:ldof
    for j=1:ldof
        KE(i,j) = problem.B(@(x, d)problem.basis_fun(x, i, d), @(x, d)...
            problem.basis_fun(x, j, d), @(x, d)evaluateActualTemperature(x, problem, actualTemperature, d), X1, X2);
    end
end

end

