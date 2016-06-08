function [ reducedBasisFuns, PODProblem ] = reducedBasis( x, iMode, derivative, coarseProblem, enrichedElementCoords )
%% Reduce Order Basis 
%generate the reduced basis used to enrich the solution space of the global
%coarse mesh. 

funs = {@getReducedBasis, @getReducedBasisDerivatives};

PODProblem = generateLocalInterpolationProblem(mapGlobalToLocal(enrichedElementCoords, enrichedElementCoords(1), enrichedElementCoords(end)));
reducedBasisFuns= funs{derivative+1}(x, iMode, coarseProblem, PODProblem );

end

function F = getReducedBasis( x, iMode, problem, PODProblem )

%Evaluate the POD basis function F(x). 

PODBasis = getPODBasis( x, problem, PODProblem, 0, iMode );
F =  PODBasis;

end

function gradF = getReducedBasisDerivatives( x, iMode, problem, PODProblem )

%Evaluate the POD basis function derivative gradF(x)...
PODBasisDerivative = getPODBasis( x, problem, PODProblem, 1, iMode );
gradF = PODBasisDerivative;

end

function PODBasis  = getPODBasis( x, problem, PODProblem, derivative, iMode )    

    PODBasis = zeros(size(x));
    X1 = PODProblem.coords(1);
    X2 = PODProblem.coords(2);

    PODBasis(x>=X1 & x<=X2) = interpolateElements(x(x>=X1 & x<=X2), problem, 1, PODProblem, derivative, iMode);
    
    for e=2:PODProblem.N
        X1 = PODProblem.coords(e);
        X2 = PODProblem.coords(e+1);
        
        PODBasis(x>X1 & x<=X2) = interpolateElements(x(x>X1 & x<=X2), problem, e, PODProblem, derivative, iMode);
    end
end

function r = interpolateElements(x, problem, element, PODProblem, derivative, iMode)  

    X1 = PODProblem.coords(element);
    X2 = PODProblem.coords(element+1);
    ldof = 2;

    r = zeros(size(x));
    for i=1:ldof
        r=r+problem.reductionOperator(PODProblem.LM(element,i), iMode).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, derivative);
    end
    
%     if derivative == 0
%         r = r.* sqrt(PODProblem.F_map(X1, X2));
%     else
%         r = r.* sqrt(PODProblem.B_map(X1, X2));
%     end

    r = r .* (2/(X2-X1)) ^ derivative;
    
end


function interpolationProblemStruct = generateLocalInterpolationProblem(coords)
    
    %Number of sub-elements
    nDOFS = size(coords, 2);
    N = nDOFS - 1;
    
    %Location Map
    PODCoefficientsLM = zeros(nDOFS, 2);
    
    %Forward and Backward map
    B_map = @(X1, X2) 2/(X2 - X1);
    F_map = @(X1, X2) (X2 - X1)/2;
    
    for i=1:nDOFS
       PODCoefficientsLM(i,1) = i;
       PODCoefficientsLM(i,2) = i+1;
    end

    interpolationProblemStruct = struct('N', N, 'nDOFS', nDOFS, 'LM', PODCoefficientsLM, 'coords', coords,...
        'B_map', B_map, 'F_map', F_map);
end
