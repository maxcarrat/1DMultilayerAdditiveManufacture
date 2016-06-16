function [ eXtendedBasisFuns ] = xFEMBasis( x, i, iMode, derivativeXBasis, derivativeLinearBasis, coarseProblem, enrichedElementCoords )
%% X-FEM Basis 
%generate the X-FEM basis used to enrich the solution space of the global
%coarse mesh. The reduced basis obtained via POD are interpolated using the
%PU basis functions of the coares/global problem, such that the enrichment
%is kept local.

if derivativeXBasis==0 && derivativeLinearBasis==0
    eXtendedBasisFuns= getXtendedBasis(x, i, iMode, coarseProblem, enrichedElementCoords );
elseif derivativeXBasis==0 && derivativeLinearBasis==1
    eXtendedBasisFuns= getXtendedBasisLinearDerivative(x, i, iMode, coarseProblem, enrichedElementCoords );
elseif derivativeXBasis==1 && derivativeLinearBasis==0
    eXtendedBasisFuns= getXtendedBasisXBasisDerivative(x, i, iMode, coarseProblem, enrichedElementCoords );
else
    eXtendedBasisFuns= getXtendedBasisTotalDerivatives(x, i, iMode, coarseProblem, enrichedElementCoords );
end
       
end

function Phi = getXtendedBasis( x, i, iMode, problem, enrichedElementCoords )

%Evaluate the enriched basis function Phi(x) using the approach presented in Fries
%and Belytschko "The extended/generalized finite element method: 
%An overview of the method and its applications"(2010). 

Xi1 = enrichedElementCoords(1);
Xi2 = enrichedElementCoords(end);

%Evaluate the ithMode basis function F(x)
if(i==1)
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    F = problem.rbBasis_fun( x, iMode, 0.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 0.0, problem, enrichedElementCoords );
else
    nodeLocalCoord = mapGlobalToLocal(Xi2,Xi1,Xi2);
    F = problem.rbBasis_fun( x, iMode, 0.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 0.0, problem, enrichedElementCoords );
end

%Evaluate the ith PU basis function N(x)
N = problem.basis_fun(x, i, 0.0);

Phi =  N .* F; %* 0.5*(enrichedElementCoords(end)-enrichedElementCoords(end-1)); 

end

function PhiDerivative = getXtendedBasisLinearDerivative( x, i, iMode, problem, enrichedElementCoords )

Xi1 = enrichedElementCoords(1);
Xi2 = enrichedElementCoords(end);

%Evaluate the ithMode basis function F(x)
if(i==1)
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    F = problem.rbBasis_fun( x, iMode, 0.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 0.0, problem, enrichedElementCoords );
else
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    F = problem.rbBasis_fun( x, iMode, 0.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 0.0, problem, enrichedElementCoords );
end

%Evaluate the ith PU basis function derivative B(x)
B = problem.basis_fun(x, i, 1.0);

PhiDerivative =  B .* F; %* 0.5*(enrichedElementCoords(end)-enrichedElementCoords(end-1)); 

end

function PhiDerivative = getXtendedBasisXBasisDerivative( x, i, iMode, problem, enrichedElementCoords )

Xi1 = enrichedElementCoords(1);
Xi2 = enrichedElementCoords(end);

%Evaluate the ithMode basis function derivative gradF(x)
if(i==1)
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    gradF = problem.rbBasis_fun( x, iMode, 1.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 1.0, problem, enrichedElementCoords );
else
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    gradF = problem.rbBasis_fun( x, iMode, 1.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 1.0, problem, enrichedElementCoords );
end

%Evaluate the ith PU basis function N(x)
N = problem.basis_fun(x, i, 0.0);

PhiDerivative =  N .* gradF; % * 0.5*(enrichedElementCoords(end)-enrichedElementCoords(end-1)); 

end

function PhiDerivative = getXtendedBasisTotalDerivatives(  x, i, iMode, problem, enrichedElementCoords )

Xi1 = enrichedElementCoords(1);
Xi2 = enrichedElementCoords(end);

%Evaluate the ithMode basis function derivative gradF(x)
if(i==1)
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    gradF = problem.rbBasis_fun( x, iMode, 1.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 1.0, problem, enrichedElementCoords );
else
    nodeLocalCoord = mapGlobalToLocal(Xi1,Xi1,Xi2);
    gradF = problem.rbBasis_fun( x, iMode, 1.0, problem, enrichedElementCoords ) -...
        problem.rbBasis_fun( nodeLocalCoord, iMode, 1.0, problem, enrichedElementCoords );
end

%Evaluate the ith PU basis function B(x)
B = problem.basis_fun(x, i, 1.0);

PhiDerivative =  B .* gradF;

end
