function globalCoordinates = mapParametricToGlobal( parametricCoordinate, problem )
%MAPPARAMETRICTOGLOBAL Map the parametric coordinates onto the global space

% N = zeros(length(parametricCoordinate), length(problem.knotVector)-problem.p-1);
% 
% for k = 1:length(parametricCoordinate)
%     [N(k,:),~] = BsplinesShapeFunctionsAndDerivatives( parametricCoordinate(k), problem.p, problem.knotVector );
% end
% 
% globalCoordinates = problem.coords * N';

globalCoordinates = zeros(length(parametricCoordinate), 1);

for k=1:length(parametricCoordinate)
    globalCoordinates(k) = parametricCoordinate(k) * (problem.coords(end)-problem.coords(1)) + problem.coords(1);
end

end

