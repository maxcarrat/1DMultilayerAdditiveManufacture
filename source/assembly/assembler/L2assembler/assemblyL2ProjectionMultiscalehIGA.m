function [ M, f ] = assemblyL2ProjectionMultiscalehIGA( solutionCoefficients, baseProblem,...
    overlayProblem, previousBaseProblem, integrationOrder, initialTemperature )
%ASSEMBLYL2PROJECTIONMULTISCALEPODXIGA assemble the matrix M and the vector f for
%the L2 projection system
%Input:
%solutionCoefficients = coefficients of the solutions
%baseProblem = poisson problem struct of the base mesh
%overlayProblem = poisson problem struct of the overlay mesh
%previousBaseProblem = problem of the previous layer base mesh
%layerLength = length of the layer
%initialtemperature = initial temperature of the powder
%Output:
%M = L2 projection matrix
%f = L2 projection rhs vector

%% Allocate matrices
%base  matrix
M_bb = zeros(baseProblem.gdof, baseProblem.gdof);
%overlay matrix
M_oo = zeros(overProblem.gdof, overlayProblem.gdof);
%coupling matrix
M_ob = zeros(baseProblem.gdof, overlayProblem.gdof);
%base  vector
f_b = zeros(baseProblem.gdof, 1);
%overlay vector
f_o = zeros(overlayProblem.gdof, 1);

%gauss points and weights
[rGP, wGP] = gaussPoints( integrationOrder );

numberOfIntegrationPoints = length(rGP);
integrationDomain = linspace(-1, 1, overlayProblem.gdof);
subDomainShapeFunctionCoefficients = linspace(0, 1, overlayProblem.gdof);


%loop over all base elements
for e=1:baseProblem.N
    
    %local dofs of the element
    ldof = baseProblem.p + 1;
    
    %knot span left and right end
    Xp1 = baseProblem.knotVector(e + baseProblem.p);
    Xp2 = baseProblem.knotVector(e + baseProblem.p + 1);
    
  
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            %map GPs onto parametric space
            localCoords = mapParentToLocal(rGP(iGP), Xp1, Xp2);
            
            %map GPs onto global coordinates system
            globalCoords = mapParentToGlobal(rGP(iGP), Xp1, Xp2, baseProblem, e);
            
            %evaluate element basis functions and derivatives
            [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords, baseProblem.p, baseProblem.knotVector);
            
            %Jacobian from parametric to global space
            JacobianX_Xi = B(baseProblem.LM(e, 1:ldof)) * baseProblem.coords(baseProblem.LM(e, 1:ldof))';
            %determinant of the Jacobian
            detJacobianX_Xi = norm(JacobianX_Xi);
            
            %Evaluate temperature@GP
            %if last element of the mesh in previous time step
            if globalCoords > baseProblem.knotVector(end - baseProblem.p - 2)
                for integrationSubDomainIndex =...
                        1 : overlayProblem.gdof-1
                    
                    if rGP(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                            rGP(iGP) > integrationDomain(integrationSubDomainIndex)
                        temperatureAtGaussPoint = evaluateTemperatureSubElements(integrationSubDomainIndex,...
                            rGP(iGP), baseProblem, previousBaseProblem, solutionCoefficients,...
                            subDomainShapeFunctionCoefficients, e, integrationDomain);
                    end
                end
            else
                temperatureAtGaussPoint = evaluateTemperature(e,localCoords, baseProblem,...
                    solutionCoefficients);
            end
            %% Integrate IGA block
            %rhs
            f_b(baseProblem.LM(e,1:ldof)) = f_b(baseProblem.LM(e,1:ldof)) + N(baseProblem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                * wGP(iGP) * detJacobianX_Xi * baseProblem.F_map(Xp1, Xp2);
            
            %projection matrix
            M_bb(baseProblem.LM(e, 1:ldof), baseProblem.LM(e, 1:ldof)) = M_bb(baseProblem.LM(e, 1:ldof), baseProblem.LM(e, 1:ldof)) +...
                (N(baseProblem.LM(e, 1:ldof))' * N(baseProblem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * baseProblem.F_map(Xp1, Xp2);
            
        end
end

%integration domain for BSpline to be integarted onto the overlay domain
integrationDomain = linspace(-1, 1, overlayProblem.N + 1);
subDomainShapeFunctionCoefficients = linspace(0, 1, overlayProblem.N + 1);

%loop over all overlay elements
for e=1:overlayProblem.N
    
    %local dofs of the element
    ldof = overlayProblem.p + 1;
    
    %nodal coordinates
    Xn1 = overlayProblem.coords(e);
    Xn2 = overlayProblem.coords(e + 1);
        
    %knot span left and right end
    Xp1 = baseProblem.knotVector(end-1-baseProblem.p);
    Xp2 = baseProblem.knotVector(end - baseProblem.p);
    
  
        % Gauss integration
        for iGP = 1:numberOfIntegrationPoints
            
            %evaluate element basis functions and derivatives
            [N, B] = BSplinesShapeFunctionsAndDerivativesSubElements( rGP(iGP), e, subDomainShapeFunctionCoefficients,...
                Xp1, Xp2, integrationDomain(e), integrationDomain(e+1), baseProblem );
            
            %Jacobian from parametric to global space
            JacobianX_Xi = B(baseProblem.LM(end, 1:ldof)) * baseProblem.coords(baseProblem.LM(end, 1:ldof))';
            %determinant of the Jacobian
            detJacobianX_Xi = norm(JacobianX_Xi);
            
            %Evaluate temperature@GP
            temperatureAtGaussPoint = initialTemperature;
            
            %% Integrate FEM block
            %rhs
            f_o(overlayProblem.LM(e,1:ldof)) = f_o(overlayProblem.LM(e,1:ldof)) + N(overlayProblem.LM(e, 1:ldof))' * temperatureAtGaussPoint...
                * wGP(iGP) * overlayProblem.F_map(Xp1, Xp2);
            
            %projection matrix
            M_oo(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) = M_oo(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) +...
                (N(overlayProblem.LM(e, 1:ldof))' * N(overlayProblem.LM(e, 1:ldof))) * wGP(iGP) * overlayProblem.F_map(Xn1, Xn2);
            
            %Integrate Coupling matrix
            M_ob(overlayProblem.LM(e, 1:ldof), overlayProblem.LCoupling(e, 1:ldof)) = M_ob(overlayProblem.LCoupling(e, 1:ldof), overlayProblem.LCoupling(e, 1:ldof)) +...
                (N(overlayProblem.LM(e, 1:ldof))' * N(overlayProblem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * overlayProblem.F_map(Xp1, Xp2);            
            
        end
end

M = [M_bb, M_ob'; M_ob, M_oo];
f = [f_b, f_o];

end


function [ projectedCoefficients ] = evaluateTemperature(e, xLocal,...
    problem, solutionCoefficients)
%%EVALUATETEMPERATURE project the previous solution onto the element.
%Input:
%e = element index
%xGlobal = global coordinates of the point to be evaluated
%xLocal = parametric coordinates of the point to be evaluated
%problem = poisson problem struct
%layerLength = length of the layer
%initialtemperature = initial temperature of the powder
%Output:
%projectedCoefficients = temperature values 

%% Initialize variables
numberOfProjectionPoints = length(xLocal);
projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
N = zeros(length(xLocal), 2);
localCoords = xLocal;

%loop over points and evaluate the BSplines at that point
for k=1:length(xLocal)
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(localCoords(k),problem.p, problem.knotVector);
end

%assembly projection operator
for i=1:length(xLocal)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end

projectedCoefficients = projectionOperator(:, problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));

end

function [ projectedCoefficients ] = evaluateTemperatureSubElements(integrationDomainIndex,...
    x, baseProblem, previousBaseProblem, solutionBaseCoefficients, solutionOverlayCoefficients, ...
    subDomainShapeFunctionCoefficients, element, integrationDomain)
% EVALUATETEMPERATURESUBELEMENTS project the previous solution onto the element.
%Input:
%integrationDomainIndex = index of the integration domain element
%x = projection local coordinates
%baseProblem = base mesh problem struct
%previousBaseProblem = previous layer mesh problem struct
%solutionBaseCoefficients = coefficients of the base mesh solutions
%solutionOverlayCoefficients = coefficients of the overlay mesh solutions
%subDomainShapeFunctionCoefficients = coefficients of the base mesh shape 
%functions (used only if p=1 !!!)
%element = index of the element
%integrationDomain = subDOmain for the integration of the BSplines basis
%Output:
%projectedCoefficients = temperature values 

numberOfProjectionPoints = length(x);

%projection operator has :
% #rows = number of projection points
% #columns = base elemnts order +1 + overlay elements order + 1
projectionOperator = zeros(numberOfProjectionPoints, baseProblem.p + 1 + 2);

solutionCoefficients = [solutionBaseCoefficients; solutionOverlayCoefficients];

N = zeros(length(x), problem.p + 1);
F = zeros(length(x), 2);

XiParametric1 = problem.knotVector( element + problem.p );
XiParametric2 = problem.knotVector( element + 1 + problem.p );

Xi1 = integrationDomain( integrationDomainIndex );
Xi2 = integrationDomain( integrationDomainIndex + 1 );

for k=1:length(x)    
    [N(k,:), ~] = BSplineShapeFunctionsAndDerivativesComposedIntegration(x(k), integrationDomainIndex,...
        subDomainShapeFunctionCoefficients, XiParametric1, XiParametric2, Xi1, Xi2, previousBaseProblem);
    [F(k,:), ~] = shapeFunctionsAndDerivatives(x(k));
end

for i=1:length(x)
    projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,1:end), F(i,:)];
end

projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * solutionCoefficients(problem.N-1:end);
end

