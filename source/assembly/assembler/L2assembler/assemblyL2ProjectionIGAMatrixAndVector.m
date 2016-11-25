function [ M, f ] = assemblyL2ProjectionIGAMatrixAndVector( solutionCoefficients, problem, integrationOrder,...
     layerLength, initialTemperature )
%ASSEMBLYL2PROJECTIONMATRIXANDVECTOR assemble the matrix M and the vector f for
%the L2 projection system
%Input:
%solutionCoefficients = coefficients of the solutions
%problem = poisson problem struct
%layerLength = length of the layer
%initialtemperature = initial temperature of the powder

%% Allocate matrices
%IGA  matrix
M_IGA = zeros(problem.gdof,problem.gdof);
%IGA  vector
f_IGA = zeros(problem.gdof, 1);

%gauss points and weights
[rGP, wGP] = gaussPoints( integrationOrder );

numberOfIntegrationPoints = length(rGP);

%loop over all elements
for e=1:problem.N
    
    %local dofs of the element
    ldof = problem.p + 1;
    
    %knot span left and right end
    Xi1 = problem.knotVector(e + problem.p);
    Xi2 = problem.knotVector(e + problem.p + 1);
    
  
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            %map GPs onto parametric space
            localCoords = mapParentToLocal(rGP(iGP), Xi1, Xi2);
            
            %map GPs onto global coordinates system
            globalCoords = mapParametricToGlobal(localCoords, problem);
            
            %evaluate element basis functions and derivatives
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);
            
            %Jacobian from parametric to global space
            JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';
            
            %determinant of the Jacobian
            detJacobianX_Xi = norm(JacobianX_Xi);
            
            %Evaluate temperature@GP
            temperatureAtGaussPoint = evaluateTemperature(e, globalCoords, localCoords, problem,...
                solutionCoefficients, layerLength, initialTemperature);
            
            %% Integrate IGA block
            %external heat source
            f_IGA(problem.LM(e,1:ldof)) = f_IGA(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
            %Capacity matrix
            M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M_IGA(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
            
        end
end

M = M_IGA;
f = f_IGA;


end


function [ projectedCoefficients ] = evaluateTemperature(e, xGlobal, xLocal,...
    problem, solutionCoefficients, layerLength, initialTemperature)
% EVALUATETEMPERATURE project the previous solution onto the element.
%Input:
%e = element index
%xGlobal = global coordinates of the point to be evaluated
%xLocal = parametric coordinates of the point to be evaluated
%problem = poisson problem struct
%layerLength = length of the layer
%initialtemperature = initial temperature of the powder

%% Initialize variables
numberOfProjectionPoints = length(xLocal);
projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
N = zeros(length(xLocal), 2);
localCoords = xLocal;

%loop over points and evaluate the BSplines at that point
for k=1:length(xLocal)
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(localCoords(k),problem.p, problem.knotVector);
    projectionOperator(k,1:size(N,2)) = N(k,:);
end

%if in the new layer set to powder temperature, project the results of the
%previous mesh otherwise
if xGlobal < problem.coords(end) && xGlobal >= ( problem.coords(end) - layerLength)
    projectedCoefficients = initialTemperature;
else
    projectedCoefficients = projectionOperator(:, problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));
end

end

