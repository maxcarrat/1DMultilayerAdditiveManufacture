function [ reducedBasisFuns, localProblem ] = localIntegrationBasis( x, i, derivative, coarseProblem, enrichedElementCoords )
%% Local Integration Basis
%integrate the linear shape functions of the coarse mesh on the local
%enriched mesh

funs = {@getReducedBasis, @getReducedBasisDerivatives};

localProblem = generateLocalInterpolationProblem(mapGlobalToLocal(enrichedElementCoords, enrichedElementCoords(1), enrichedElementCoords(end)));
reducedBasisFuns= funs{derivative+1}(x, i, coarseProblem, localProblem );

end

function F = getReducedBasis( x, i, problem, PODProblem )

%Evaluate the POD basis function F(x).

PODBasis = getlocalBasis( x, i, problem, PODProblem, 0 );
F =  PODBasis;

end

function gradF = getReducedBasisDerivatives( x, i, problem, localProblem )

%Evaluate the POD basis function derivative gradF(x)...
PODBasisDerivative = getlocalBasis( x, i, problem, localProblem, 1 );
gradF = PODBasisDerivative;

end

function localBasis  = getlocalBasis( x, i, problem, localProblem, derivative )

localBasis = zeros(size(x));
X1 = localProblem.coords(1);
X2 = localProblem.coords(2);

localBasis(x>=X1 & x<=X2) = interpolateElements(x(x>=X1 & x<=X2), i, problem, 1, localProblem, derivative);

for e=2:localProblem.N
    X1 = localProblem.coords(e);
    X2 = localProblem.coords(e+1);
    
    localBasis(x>X1 & x<=X2) = interpolateElements(x(x>X1 & x<=X2), i, problem, e, localProblem, derivative);
end
end

function r = interpolateElements(x, index, problem, element, localProblem, derivative)

X1 = localProblem.coords(element);
X2 = localProblem.coords(element+1);
X = [X1, X2];
ldof = 2;

r = zeros(size(x));

for i=1:ldof
    r=r+ problem.basis_fun(X(i), index, 0.0) * problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, derivative);
end


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
