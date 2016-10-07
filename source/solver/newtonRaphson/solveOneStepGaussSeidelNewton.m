function [ baseSolution, overlaySolution, convergenceFlag ] = solveOneStepGaussSeidelNewton( baseProblem,...
    overlayProblem, time, timeStepSize, integrationOrder, overlayIntegrationOrder, tolerance,...
    maxNumberOfIterations, oldBaseSolution, oldOverlaySolution, eXtendedFlag)
%SOLVEONESTEPGAUSSSEIDELNEWTON returns the temperature solution of the base
%and the overlay mesh. The nonlinear system is solved using the one-step
%Gauss-Seidel-Newton scheme, which allows to solve a multiscale hp-d
%nonlinear problem with a single newton iteration and a Gauss-Seidel loop
%Input:
%baseProblem = IBV problem struct of the base mesh
%overlayProblem = IBV problem struct of the overlay mesh
%time = actual time of the process
%timeStepSize = size of the temporal discretization
%integrationOrder = number of quadrature point in the base mesh
%overlayIntegrationOrder = number of qudrature points in the overlay mesh
%tolerance = convergence tolerance 
%maxNumberOfIterations = max number of Gauss-Seidel iterations
%oldBaseSolution = last converged solution of the base mesh
%oldOverlaySolution = last converged solution of the overlay mesh
%eXtendedFlag = true if teh mesh is eXtended , false otherwise
%Output:
%baseSolution = solution coefficients of the base mesh
%overlaySolutions = solutions coefficients of the overlay mesh
%convergenceFlag = check if the solution converged or not

%base mesh
baseSolution = oldBaseSolution;
lastConvergedBaseSolution = oldBaseSolution;

%overlay mesh
overlaySolution = oldOverlaySolution;
lastConvergedOverlaySolution = oldOverlaySolution;

convergenceFlag = 1;

%Gauss-Seidel loop
for i=1:maxNumberOfIterations
    
   %% Compute base mesh solution increment --------------------------------
   %assembly base mesh ...
   [M_bb, Mprime_bb, K_bb, Kprime_bb, f_b] = assemblyMultiscaleBaseIGA(baseProblem, overlayProblem, time, integrationOrder,...
       baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
   
   if strcmp(eXtendedFlag, 'false')
       %... assembly coupling terms...
       [M_ob, K_ob] = assemblyMultiscaleCouplingFEM(baseProblem, overlayProblem, time, intrgationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
   else
       %... assembly coupling terms...
       [M_ob, K_ob] = assemblyMultiscaleCouplingPODXFEM(baseProblem, overlayProblem, time, intrgationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
   end
   %... and apply boundary comditions.
   [K_bb, Kprime_bb, f_b] = applyGlobalBCs(baseProblem, K_bb, Kprime_bb, f_b);
   
   %residuum
   r = f_b * timeStepSize - (K_bb * (baseSolution) * timeStepSize +...
       K_ob' * (overlaySolution) * timeStepSize ) +...
       M_bb * (lastConvergedBaseSolution) - M_bb * (baseSolution) -...
       M_ob' * (lastConvergedOverlaySolution) - M_ob' * (overlaySolution);
   
   %jacobian of the residuum
   J = ( K_bb + Kprime_bb ) * timeStepSize +...
       (M_bb + Mprime_bb );
   
   %solve
   baseSolutionIncrement = J\r;
   
   %update base solution
   baseSolution = baseSolution + baseSolutionIncrement;
   
   %% Compute overlay mesh solution increment --------------------------------
  %if the mesh is enriched:
   if strcmp(eXtendedFlag, 'false')
       %assembly overlay mesh ...
       [M_oo, Mprime_oo, K_oo, Kprime_oo, f_o] = assemblyMultiscaleOverlayFEM(baseProblem, overlayProblem, time, overlayIntegrationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
       
       %... assembly coupling terms ...
       [M_ob, K_ob] = assemblyMultiscaleCouplingFEM(baseProblem, overlayProblem, time, overlayIntegrationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
   else
   % ... otherwise:
       %assembly overlay mesh ...
       [M_oo, Mprime_oo, K_oo, Kprime_oo, f_o] = assemblyMultiscaleOverlayPODXFEM(baseProblem, overlayProblem, time, overlayIntegrationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
       
       %... assembly coupling terms...
       [M_ob, K_ob] = assemblyMultiscaleCouplingPODXFEM(baseProblem, overlayProblem, time, overlayIntegrationOrder,...
           baseSolution, overlaySolution, lastConvergedBaseSolution, lastConvergentOverlaySolution);
   end
   %... and apply boundary comditions.
   [K_oo, Kprime_oo, f_o] = applyGlobalBCs(baseProblem, K_oo, Kprime_oo, f_o);
   
   %residuum
   r = f_o * timeStepSize - (K_oo * (baseSolution) * timeStepSize +...
       K_ob' * (overlaySolution) * timeStepSize ) +...
       M_oo * (lastConvergedBaseSolution) - M_oo * (baseSolution) -...
       M_ob' * (lastConvergedOverlaySolution) - M_ob' * (overlaySolution);
   
   %jacobian of the residuum
   J = ( K_oo + Kprime_oo ) * timeStepSize +...
       (M_oo + Mprime_oo );
   
   %solve
   overlaySolutionIncrement = J\r;

   %update base solution
   overlaySolution = overlaySolution + overlaySolutionIncrement;
   
   %% Check convergence in L2-norm ----------------------------------------
   %total solution
   solution = baseSolution + overlaySolution;
   resNorm = sqrt(solution'*solution);
   
   %check convergence
   if resNorm <= tolerance
       formatSpec = '\t\t Solution converge after %1.0f iterations\n';
       fprintf(formatSpec, i);
       convergenceFlag = 0;
       return
   end
   
   formatSpec = '\t\t Residuum Norm %1.5f \n';
   fprintf(formatSpec, resNorm);
       
end

formatSpec = '\t\t Warning: Solution did not converge after %1.0f iterations!!!\n';
fprintf(formatSpec, maxNumberOfIterations);
   
end


