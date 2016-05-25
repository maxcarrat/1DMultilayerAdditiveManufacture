function [ reducedBasisFuns ] = reducedBasis( x, i, iMode, derivative, coarseProblem, enrichedElementCoords )
%% Reduce Order Basis 
%generate the reduced basis used to enrich the solution space of the global
%coarse mesh

funs = {@getReducedBasis, @getReducedBasisDerivatives};
reducedBasisFuns= funs{derivative+1}(x, i, iMode, coarseProblem, enrichedElementCoords );

end

function integratedFunction = getReducedBasis( x, i, iMode, problem, enrichedElementCoords )
integratedFunction = zeros(1, size(x,2));

for k=1:size(x)
    for j=1:size(enrichedElementCoords,2)-1
        
        if(mapLocalToGlobal(x(k), enrichedElementCoords(1), enrichedElementCoords(end))< enrichedElementCoords(j+1)) &&...
                (mapLocalToGlobal(x(k), enrichedElementCoords(1), enrichedElementCoords(end))>= enrichedElementCoords(j))
            
            integratedFunction = problem.basis_fun(x, i, 0)*problem.reductionOperator(j, iMode);
        
        end
    end
end

end

function integratedFunction = getReducedBasisDerivatives( x, i, iMode, problem, enrichedElementCoords )

integratedFunction = zeros(1, size(x,2));

for k=1:size(x)
    for j=1:size(enrichedElementCoords,2)-1
        
        if(mapLocalToGlobal(x(k), enrichedElementCoords(1), enrichedElementCoords(end))< enrichedElementCoords(j+1)) &&...
                (mapLocalToGlobal(x(k), enrichedElementCoords(1), enrichedElementCoords(end))>= enrichedElementCoords(j))
            
            integratedFunction = problem.basis_fun(x, i, 0)*problem.reductionOperator(j, iMode);
            
        end
    end
end
end
