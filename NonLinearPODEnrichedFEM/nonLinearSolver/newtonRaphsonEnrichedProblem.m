function [ solution ] = newtonRaphsonEnrichedProblem( oldSolution, problem, numberOfModesSupports,...
    tolerance, maxNumberOfIterations, timeStepSize )
%NEWTONRAPHSONENRICHEDPROBLEM Evaluate the solution cofficients at each time step using
%Newton-Raphson algorithm of the POD-Enriched problem
%   oldSolution = solution coefficients of the previous iteration
%   problem = Poisson transient non-linear problem
%   tolerance = convergence tolerance
%   maxNumberOfIteration = max number of iteration of the algorithm
%   timeStepSize = size of the time increment at each time step

% solution = oldSolution;
% solutionRB = oldSolution(end-numberOfModesSupports:end);
% lastConvergedSolution = oldSolution(end-numberOfModesSupports:end);

solution = oldSolution;
solutionRB = zeros(numel(oldSolution), 1);
lastConvergedSolution = oldSolution;

for i=1:maxNumberOfIterations
    
    %% Assembly and apply BCs
    
    %Assembly the local reduced basis
    [M, K, f] = assemblyLocalNonLinearProblem(problem, solution);
    
    %Apply BCs
    constrainedTemperature = problem.dirichlet_bc(2, 2);
    constrainedNode = problem.coords(end);
%     [K, f] = applyWeakDirichletBCs(problem, constrainedNode, constrainedTemperature, K, f);
    K(2,2) = K(2,2) + problem.penalty;
    f(2) = f(2) + constrainedTemperature*problem.penalty;

    %% Solve
    %residuum
    r = f * timeStepSize - K * (solutionRB) * timeStepSize +...
        M * (lastConvergedSolution) - M * (solutionRB);
    
    %jacobian of the residuum
    J = K * timeStepSize + M;
    
    solutionIncrement = J\r;
    resNorm = sqrt(abs(r'*solutionIncrement));
    
    %% Check convergence
    if resNorm <= tolerance
        formatSpec = '\t\t Solution converge after %1.0f iterations\n';
        fprintf(formatSpec, i);
        solutionRB = solutionRB + solutionIncrement;
        solution = solutionRB;
        return
    end
    
    formatSpec = '\t\t Residuum Norm %1.5f \n';
    fprintf(formatSpec, resNorm);
    
    solutionRB = solutionRB + solutionIncrement;    
    solution = solutionRB;

end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);


end