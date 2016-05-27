function [ solution ] = newtonRaphsonGlobalProblem( oldSolution, problem, tolerance, maxNumberOfIterations, timeStepSize)
%NEWTONRAPHSONGLOBALPROBLEM Evaluate the solution cofficients at each time step using
%Newton-Raphson algorithm
%   oldSolution = solution coefficients of the previous iteration
%   problem = Poisson transient non-linear problem
%   tolerance = convergence tolerance
%   maxNumberOfIteration = max number of iteration of the algorithm
%   timeStepSize = size of the time increment at each time step

solution = oldSolution;
lastConvergedSolution = oldSolution;

for i=1:maxNumberOfIterations
    
   %assembly and apply BCs
   [M, K, f] = assemblyNonLinearProblem(problem, solution);
   [M, K, f] = applyBCsNonLinearProblem(M, K, f, problem);
   
   %residuum
   r = f * timeStepSize - K * (solution) * timeStepSize +...
       M * (lastConvergedSolution) - M * (solution);
   
   %jacobian of the residuum
   J = K * timeStepSize + M;
   
   solutionIncrement = J\r;
   resNorm = norm(solutionIncrement);
   
   %check convergence
   if resNorm <= tolerance
       formatSpec = '\t\t Solution converge after %1.0f iterations\n';
       fprintf(formatSpec, i);
       solution = solution + solutionIncrement;
       return
   end
   
   formatSpec = '\t\t Residuum Norm %1.5f \n';
   fprintf(formatSpec, resNorm);
   
   solution = solution + solutionIncrement;
    
end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);
   
end

