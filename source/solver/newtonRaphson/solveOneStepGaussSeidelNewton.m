function [ baseSolution, overlaySolution, convergenceFlag ] = solveOneStepGaussSeidelNewton( initialOverlayProblem,...
    overlayProblem, baseProblem, time, timeStepSize, integrationOrder, overlayIntegrationOrder, tolerance,...
    maxNumberOfIterations, oldBaseSolution, oldOverlaySolution, initialTemperature, eXtendedFlag)
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
lastConvergentBaseSolution = oldBaseSolution;

%overlay mesh
% overlaySolution = zeros(size(oldOverlaySolution,1),size(oldOverlaySolution,2));
% lastConvergentOverlaySolution = zeros(size(oldOverlaySolution,1),size(oldOverlaySolution,2));

overlaySolution = oldOverlaySolution;
lastConvergentOverlaySolution = oldOverlaySolution;

convergenceFlag = 1;

%Gauss-Seidel loop
for i=1:maxNumberOfIterations
    
    %% Compute base mesh solution increment --------------------------------
    %assembly base mesh ...
    
    if strcmp(eXtendedFlag, 'false')
        %... assembly linear system...
        [M_bb, Mprime_bb, K_bb, Kprime_bb, f_b] = assemblyMultiscaleBaseIGA( overlayProblem, baseProblem, time, integrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        [M_bo, K_bo] = assemblyMultiscaleCouplingBaseFEM( overlayProblem, baseProblem, time, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
    else
        %... assembly linear system...
        [M_bb, Mprime_bb, K_bb, Kprime_bb, f_b] = assemblyMultiscaleBaseXIGA( overlayProblem, baseProblem, time, integrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution, initialTemperature);
        [M_bo, K_bo] = assemblyMultiscaleCouplingPODXFEM( initialOverlayProblem, overlayProblem, baseProblem, time, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        K_bo = K_bo';
        M_bo = M_bo';
    end
    %... and apply boundary comditions.
    [K_bb, Kprime_bb, K_bo, M_bo, f_b] = applyBaseMultiscaleBCs(baseProblem, K_bo, K_bb, Kprime_bb, M_bo, f_b);
    
    %residuum
    deltaBaseSolution = (lastConvergentBaseSolution - baseSolution);
    deltaOverlaySolution =  -(lastConvergentOverlaySolution - overlaySolution);
    
    r = f_b * timeStepSize - K_bb * (baseSolution) * timeStepSize -...
        - M_bb * (deltaBaseSolution) - K_bo * (overlaySolution) * timeStepSize...
        - M_bo * (deltaOverlaySolution);
    
    %jacobian of the residuum
    J = ( K_bb + Kprime_bb ) * timeStepSize +...
        ( M_bb + Mprime_bb );
    
    %solve
    baseSolutionIncrement = J\r;
    
    %update base solution
    baseSolution = baseSolution + baseSolutionIncrement;
    
    %% Compute overlay mesh solution increment --------------------------------
    %if the mesh is enriched:
    if strcmp(eXtendedFlag, 'false')
        %assembly overlay mesh ...
        [M_oo, Mprime_oo, K_oo, Kprime_oo, f_o] = assemblyMultiscaleOverlayFEM( overlayProblem, baseProblem, time, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        
        %... assembly coupling terms ...
        [M_ob, K_ob] = assemblyMultiscaleCouplingOverlayFEM( overlayProblem, baseProblem, time, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        
        %... and apply boundary comditions.
        [K_oo, Kprime_oo, K_ob, M_oo, M_ob, Mprime_oo, f_o] = applyOverlayMultiscaleBCs(overlayProblem, K_oo, Kprime_oo, K_ob, M_oo, M_ob, Mprime_oo, f_o);
        
        %residuum
        deltaBaseSolution = -(lastConvergentBaseSolution - baseSolution);
        deltaOverlaySolution = (lastConvergentOverlaySolution - overlaySolution);
        
        r = f_o * timeStepSize - K_oo * (overlaySolution) * timeStepSize -...
            K_ob * (baseSolution) * timeStepSize -...
            M_oo * (deltaOverlaySolution) - M_ob * (deltaBaseSolution);
        
        %jacobian of the residuum
        J = ( K_oo + Kprime_oo ) * timeStepSize +...
            ( M_oo + Mprime_oo );
        
        %solve
%         J(1,:) = [];
%         J(end,:) = [];
%         J(:,1) = [];
%         J(:,end) = [];
%         r(1) = [];
%         r(end) = [];
        overlaySolutionIncrement = J\r;
        
        %apply Dirichlet Boundaries strongly on the solution
%         overlaySolutionIncrement = [0.0;overlaySolutionIncrement;0.0];

%         overlaySolutionIncrement(1) = 0.0;
%         overlaySolutionIncrement(end) = 0.0;
    else
        % ... otherwise:
        %assembly overlay mesh ...
        [M_oo, Mprime_oo, K_oo, Kprime_oo, f_o] = assemblyMultiscaleOverlayPODXFEM( initialOverlayProblem, overlayProblem, baseProblem, time,...
            integrationOrder, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        
        %... assembly coupling terms...
        [M_ob, K_ob] = assemblyMultiscaleCouplingPODXFEM( initialOverlayProblem, overlayProblem, baseProblem,  time, overlayIntegrationOrder,...
            overlaySolution, baseSolution, lastConvergentBaseSolution, lastConvergentOverlaySolution);
        
        %... and apply boundary comditions.
        [K_oo, Kprime_oo, K_ob, M_oo, M_ob, Mprime_oo, f_o] = applyOverlayXtendedMultiscaleBCs(overlayProblem, K_oo, Kprime_oo, K_ob, M_oo, M_ob, Mprime_oo, f_o);
        
        %residuum
        deltaBaseSolution = -(lastConvergentBaseSolution - baseSolution);
        deltaOverlaySolution = +(lastConvergentOverlaySolution - overlaySolution);
        r = f_o * timeStepSize - K_oo * (overlaySolution) * timeStepSize -...
            K_ob * (baseSolution) * timeStepSize -...
            M_oo * (deltaOverlaySolution) - M_ob * (deltaBaseSolution);
        
        %jacobian of the residuum
        J = ( K_oo + Kprime_oo ) * timeStepSize +...
            ( M_oo + Mprime_oo );
        
        %solve
%         J(1:2,:) = [];
%         J(:,1:2) = [];
%         r(1:2) = [];
        overlaySolutionIncrement = J\r;
        
        
        %apply Dirichlet Boundaries strongly on the solution
%         overlaySolutionIncrement = [0.0;0.0;overlaySolutionIncrement];

%         overlaySolutionIncrement(1) = 0.0;
%         overlaySolutionIncrement(2) = 0.0;
    end

    %update base solution
    overlaySolution = overlaySolution + overlaySolutionIncrement;
    
    %% Check convergence in L2-norm ----------------------------------------
    %total solution
    solutionIncrement = [baseSolutionIncrement; overlaySolutionIncrement];
    resNorm = sqrt(solutionIncrement'*solutionIncrement);
    
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


