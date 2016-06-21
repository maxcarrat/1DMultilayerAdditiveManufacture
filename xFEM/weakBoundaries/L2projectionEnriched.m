function [ extractedPODCoefficients ] = L2projectionEnriched(problem, previousTemperature, updatedMesh, coarseMesh,...
    modes, initialTemperature)
%L2PROJECTIONENRICHED project the previous solution onto the updated mesh at the
%new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   coarseMesh = initial/coarse mesh configuration
%   modes = number of enrichment modes
%   initialTemperature = initial temperature of the powder

ldof = 2;

projectionOperator = zeros(modes*2+2, modes*2+2);
enrichedElementCoords = linspace(updatedMesh(end-1), updatedMesh(end), 2^problem.refinementDepth+1);

[xGL, wGL] = GaussLobatto(modes*2+2);

for j=1:modes*2+2
    
    xi = xGL(j);
    wi = wGL(j);
    
    projectionOperator(j, 1) =  problem.basis_fun(xi, 1, 0.0);
    projectionOperator(j, 2) =  problem.basis_fun(xi, 2, 0.0);
    
    for i=1:ldof
        for iMode=1:modes
            projectionOperator(j, 2 + (i-1)*modes + iMode) =  problem.xFEMBasis_fun(xi, i, iMode, 0.0, 0.0, problem, enrichedElementCoords);
        end
    end
    %     projectionOperator(j, :) = wi * projectionOperator(j,:);
    
end

initialProjectedCoefficients = zeros(modes*2+2, 1);
initialProjectedCoefficients(:) = initialTemperature;
% initialProjectedCoefficients(1) = problem.dirichlet_bc(1,2);
extractedPODCoefficients = projectionOperator\initialProjectedCoefficients;


end

