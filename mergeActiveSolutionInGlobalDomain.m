function [ globalSolutionUpdated ] = mergeActiveSolutionInGlobalDomain( activeSolution, globalSolutionSize )
%MERGEACTIVESOLUTIONINGLOBALDOMAIN merge the active solutions in the global
%domain

globalSolutionUpdated = zeros(globalSolutionSize, 1);

for i=1:size(activeSolution,1)
    
    globalSolutionUpdated(i) = activeSolution(i);
    
end

end

