function [ globalTemperatureCoefficents ] = solveLocalProblem( globalTemperatureCoefficents, problem, timeStepSize, numberOfModesSupports )
%SOLVEGLOBALPROBLEM returns the nodal temperature values of the
%local/enriched problem.
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at the current time step

%% Loop over POD modes 
%The enrichment solution of each mode is added to the final enrichment
%vector

M = zeros(numberOfModesSupports+1, numberOfModesSupports+1);
K = zeros(numberOfModesSupports+1, numberOfModesSupports+1);
f = zeros(numberOfModesSupports+1, 1);

for iMode=1:problem.modes
    
    %Assembly the local reduced basis
    [M_mode, K_mode, f_mode] = assemblyLocalProblem(problem, iMode, numberOfModesSupports);
    
    M = M + M_mode;
    K = K + K_mode;
    f = f + f_mode;
end

%Apply the Dirichlet BCs at the nodes of the local enriched element

K(end, end) = K(end, end) + problem.penalty;
f(end) =  f(end) + problem.penalty * problem.dirichlet_bc(2, 2);

%Solve
previousSolution = globalTemperatureCoefficents(end-numberOfModesSupports:end);
RHS = timeStepSize * (f - K * previousSolution);
LHS = M + timeStepSize * K;
increment = LHS\RHS;

globalTemperatureCoefficents(end-numberOfModesSupports:end) = globalTemperatureCoefficents(end-numberOfModesSupports:end) + increment;

end

