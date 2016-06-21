function [ temperatureCoefficents ] = solveGlobalProblemGaussIntegration(temperatureCoefficents, problem, time, timeStepSize, integrationOrder )
%SOLVEGLOBALPROBLEM returns the nodal temperature values of the
%coarse/global problem using Backward Euler implictit scheme
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at teh current time step

%Backward Euler Scheme
[uncontrainedCapacityMatrix, unconstrainedConductivityMatrix, unconstrainedExternalSource] = assemblyLinearSystem(problem, time, integrationOrder);
[M, K, f] = applyGlobalBCs(problem, uncontrainedCapacityMatrix,...
    unconstrainedConductivityMatrix, unconstrainedExternalSource);

RHS = timeStepSize * (f - K * temperatureCoefficents);
LHS = M + timeStepSize * K;

temperatureIncrement = LHS\RHS;

temperatureCoefficents = temperatureCoefficents + temperatureIncrement;


end

