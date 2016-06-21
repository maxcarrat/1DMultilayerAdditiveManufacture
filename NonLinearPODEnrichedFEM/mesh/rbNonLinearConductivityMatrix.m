function [ K ] = rbNonLinearConductivityMatrix( problem, solution,...
    enrichedElementCoords, e)
%RBNONLINEARCONDUCTIVITYMATRIX computes the conductivity local to one element
%   problem = definition of the boundary value problem
%   solution = actual temperature value
%   iModes = index of the ith POD mode
%   enrichedElementCoords = coordinates of the refined element
%   e = element index

X1 = problem.coords(e);
X2 = problem.coords(e+1);
ldof = 2;
numberOfModes = problem.modes;

%sub-matrices
K_FEM = zeros(ldof, ldof);
K_enr = zeros(numberOfModes*ldof, numberOfModes*ldof);
K_coupling = zeros(numberOfModes*ldof, ldof);

for i=1:ldof
    for j=1:ldof
        K_FEM(i,j) = problem.B(@(x)problem.localBasis_fun(x, i, 1.0, problem, enrichedElementCoords),...
                @(x)problem.localBasis_fun(x, j, 1.0, problem, enrichedElementCoords), @(x, d)evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2);
    end
end

for i=1:ldof
    for iMode = 1:numberOfModes
        for j= 1:ldof
            for jMode = 1:numberOfModes
                K_enr((i-1)*numberOfModes + iMode,(j-1)*numberOfModes + jMode) = problem.rbB(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 1.0, problem, enrichedElementCoords ),...
                    @(x)problem.xFEMBasis_fun(x, j, jMode, 0.0, 1.0, problem, enrichedElementCoords), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2) + ...
                    problem.rbB(@(x)problem.xFEMBasis_fun(x, i, iMode, 1.0, 0.0, problem, enrichedElementCoords ),...
                    @(x)problem.xFEMBasis_fun(x, j, jMode, 0.0, 1.0, problem, enrichedElementCoords), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2) + ...
                    problem.rbB(@(x)problem.xFEMBasis_fun(x, i, iMode, 0.0, 1.0, problem, enrichedElementCoords ),...
                    @(x)problem.xFEMBasis_fun(x, j, jMode, 1.0, 0.0, problem, enrichedElementCoords), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2) + ...
                    problem.rbB(@(x)problem.xFEMBasis_fun(x, i, iMode, 1.0, 0.0, problem, enrichedElementCoords ),...
                    @(x)problem.xFEMBasis_fun(x, j, jMode, 1.0, 0.0, problem, enrichedElementCoords), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2);
            end
        end
    end
end

for i=1:ldof
    for j=1:ldof
        for jMode = 1:numberOfModes
            K_coupling((j-1)*numberOfModes + jMode, i) = problem.rbB(@(x)problem.localBasis_fun(x, i, 1.0, problem, enrichedElementCoords),...
                @(x)problem.xFEMBasis_fun(x, j, jMode, 0.0, 1.0, problem, enrichedElementCoords ), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2) + ...
                problem.rbB(@(x)problem.localBasis_fun(x, i, 1.0, problem, enrichedElementCoords),...
                @(x)problem.xFEMBasis_fun(x, j, jMode, 1.0, 0.0, problem, enrichedElementCoords ), @(x, d) evaluateLocalActualTemperature( x, problem, enrichedElementCoords, solution, d ), X1, X2);
        end
    end
end

K = [K_FEM, (K_coupling)'; K_coupling, K_enr];


end

