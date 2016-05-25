function [ enrichmentCoefficients ] = solveLocalProblem( globalTemperatureCoefficents, enrichmentCoefficients, problem, timeStepSize )
%SOLVEGLOBALPROBLEM returns the nodal temperature values of the
%local/enriched problem.
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at the current time step

%% Loop over POD modes 
%The enrichment solution of each mode is added to the final enrichment
%vector

numberOfModeSupports = 2;

for iMode=1:problem.modes
    
    %Assembly the local reduced basis
    [M, K, f] = assemblyLocalProblem(problem, iMode);
    
    %Apply the Dirichlet BCs at the nodes of the local enriched element
    
    K(1, 1) = K(1, 1) + problem.penalty;
    f(1) =  f(1) + problem.penalty * globalTemperatureCoefficents(end-1);
    
    K(end, end) = K(end, end) + problem.penalty;
    f(end) =  f(end) + problem.penalty * 0.0;
    
    %Solve
    RHS = timeStepSize * (f - K * enrichmentCoefficients);
    LHS = M + timeStepSize * K;
%     RHS = f;
%     LHS = K;
    
    enrichmentCoefficients = LHS\RHS;
    
end

%% Constrine Coefficients 
enrichmentCoefficients(1) = 0.0;
enrichmentCoefficients(end) = 0.0;

end

