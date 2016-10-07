function [ solution, convergenceFlag ] = solveMultiPhaseProblem( problem, time, timeStepSize,...
    integrationOrder, tolerance, maxNumberOfIterations, oldSolution )
%SOLVENONLINEARPROBLEMGAUSSINTEGRATION returns the nodal temperature values of the
%coarse/global problem using Backward Euler implictit scheme, if it
%converges the flag is et to 0, 1 otherwise.
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at teh current time step

solution = oldSolution;
lastConvergedSolution = oldSolution;
convergenceFlag = 1;

for i=1:maxNumberOfIterations
    
   %assembly and apply BCs
   [M, K, f] = assemblyMultiPhaseSystem(problem, time, integrationOrder,...
       solution, lastConvergedSolution);
   
   [M, K, f] = applyGlobalBCs(problem, M, K, f);
   
   %residuum
   r = f * timeStepSize - K * (solution) * timeStepSize +...
       M * (lastConvergedSolution) - M * (solution);
   
   %jacobian of the residuum
   J = K * timeStepSize + M;
   
   solutionIncrement = J\r;
   resNorm = sqrt(solutionIncrement'*solutionIncrement);
   
   %check convergence
   if resNorm <= tolerance
       formatSpec = '\t\t Solution converge after %1.0f iterations\n';
       fprintf(formatSpec, i);
       solution = solution + solutionIncrement;
       convergenceFlag = 0;
       return
   end
   
   formatSpec = '\t\t Residuum Norm %1.5f \n';
   fprintf(formatSpec, resNorm);
   
   solution = solution + solutionIncrement;
    
end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);
   
end

