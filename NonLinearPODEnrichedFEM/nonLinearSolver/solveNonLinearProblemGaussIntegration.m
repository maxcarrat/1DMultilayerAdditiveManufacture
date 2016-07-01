function [ solution ] = solveNonLinearProblemGaussIntegration( problem, time, timeStepSize, integrationOrder, tolerance, maxNumberOfIterations, oldSolution )
%SOLVENONLINEARPROBLEMGAUSSINTEGRATION returns the nodal temperature values of the
%coarse/global problem using Backward Euler implictit scheme
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at teh current time step

solution = oldSolution;
lastConvergedSolution = oldSolution;

for i=1:maxNumberOfIterations
    
   %assembly and apply BCs
   [M, K, f] = assemblyNonLinearSystem(problem, time, integrationOrder, solution);
   [M, K, f] = applyGlobalBCs(problem, M, K, f);
   
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

