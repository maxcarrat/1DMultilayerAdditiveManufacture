function [ armonicBasisFun ] = armonicBasis( x, iMode, derivative, problem, enrichedElementCoords )
%ARMONICBASIS generate the reduced basis used to enrich the solution space of the global
%coarse mesh. 

funs = {@getArmonicBasis, @getArmonicBasisDerivatives};

armonicBasisFun= funs{derivative+1}(x, iMode, problem, 0.0 );

end

function F = getArmonicBasis( x, iMode, problem, PODProblem )

%Evaluate the armonic basis function F(x). 
F =  cos(pi * (iMode) * x);

end

function gradF = getArmonicBasisDerivatives( x, iMode, problem, PODProblem )

%Evaluate the armonic basis function derivative gradF(x)...
gradF = -sin(pi * (iMode) * x);

end


