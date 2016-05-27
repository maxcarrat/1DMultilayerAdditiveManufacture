function [ solution ] = newtonRaphsonEnrichedProblem( oldSolution, problem, numberOfModesSupports,...
tolerance, maxNumberOfIterations, timeStepSize )
%NEWTONRAPHSONENRICHEDPROBLEM Evaluate the solution cofficients at each time step using
%Newton-Raphson algorithm of the POD-Enriched problem
%   oldSolution = solution coefficients of the previous iteration
%   problem = Poisson transient non-linear problem
%   tolerance = convergence tolerance
%   maxNumberOfIteration = max number of iteration of the algorithm
%   timeStepSize = size of the time increment at each time step

solution = oldSolution;
solutionRB = oldSolution(end-numberOfModesSupports:end);
lastConvergedSolution = oldSolution(end-numberOfModesSupports:end);

for i=1:maxNumberOfIterations
    
   %% Assembly and apply BCs
       % Loop over POD modes
       %The enrichment solution of each mode is added to the final enrichment
       %vector

       M = zeros(numberOfModesSupports+1, numberOfModesSupports+1);
       K = zeros(numberOfModesSupports+1, numberOfModesSupports+1);
       f = zeros(numberOfModesSupports+1, 1);

       for iMode=1:problem.modes

           %Assembly the local reduced basis
           [M_mode, K_mode, f_mode] = assemblyLocalNonLinearProblem(problem, solution, iMode, numberOfModesSupports);

           M = M + M_mode;
           K = K + K_mode;
           f = f + f_mode;
       end
       
       K(end, end) = K(end, end) + problem.penalty;
       f(end) =  f(end) + problem.penalty * problem.dirichlet_bc(2, 2);
   
   %% Solve
   %residuum
   r = f * timeStepSize - K * (solutionRB) * timeStepSize +...
       M * (lastConvergedSolution) - M * (solutionRB);
   
   %jacobian of the residuum
   J = K * timeStepSize + M;
   
   solutionIncrement = J\r;
   resNorm = norm(solutionIncrement);
   
   %% Check convergence
   if resNorm <= tolerance
       formatSpec = '\t\t Solution converge after %1.0f iterations\n';
       fprintf(formatSpec, i);
       solutionRB = solutionRB + solutionIncrement;
       solution(end-numberOfModesSupports:end) = solutionRB;
       return
   end
   
   formatSpec = '\t\t Residuum Norm %1.5f \n';
   fprintf(formatSpec, resNorm);
   
   solutionRB = solutionRB + solutionIncrement;
   solution(end-numberOfModesSupports:end) = solutionRB;
       
end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);
   
end


