function [ solution, convergenceFlag ] = solveXIGAMultiPhase( problem, time, timeStepSize, integrationOrder,...
    integrationModalOrder, tolerance, maxNumberOfIterations, oldSolution )
%XIGAMULTIPHASESOLVER returns the nodal temperature values of the
%coarse/global problem using Backward Euler implictit scheme, if it 
%converges the flag value is set to 0, 1 otherwise.
%Input:
%problem = struct which defines the initial boundary value problem
%time = actual time
%timeStepSize = size of the integrtion time step
%integrationOrder = number of quadrature points
%integrationModalOrder 0 number of quadrature points in enriched elements
%tolerance = Newton-Raphson convergence tolerance
%maxNumberOfIterations = max number of itertions in NR
%oldSolution = last convreged solution
%Output:
%solution = converged solution coefficients
%convergenceFlag = flag is 0 if the solution converged and 1 otherwise

solution = oldSolution;
lastConvergedSolution = oldSolution;
convergenceFlag = 1;

for i=1:maxNumberOfIterations
    
   %assembly and apply BCs
   [M, K_unconstrained, f_unconstrained] = assemblyMultiPhaseXIGASystem(problem, time, integrationOrder,...
       integrationModalOrder, solution, lastConvergedSolution);
   [M, K, f] = applyGlobalBCs(problem, M, K_unconstrained, f_unconstrained);
   
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
       convergenceFlag = 0;
       return
   end
   
   formatSpec = '\t\t Residuum Norm %1.5f \n';
   fprintf(formatSpec, resNorm);
   
   %update solution
   solution = solution + solutionIncrement;
    
end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);
 
end

