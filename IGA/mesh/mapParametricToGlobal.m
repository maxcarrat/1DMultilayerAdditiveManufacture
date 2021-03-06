function globalCoordinates = mapParametricToGlobal( parametricCoordinate, problem )
%MAPPARAMETRICTOGLOBAL Map the parametric coordinates onto the global space

% N = zeros(length(parametricCoordinate), length(problem.knotVector)-problem.p-1);
% 
% for k = 1:length(parametricCoordinate)
%     [N(k,:),~] = BsplinesShapeFunctionsAndDerivatives( parametricCoordinate(k), problem.p, problem.knotVector );
% end
% 
% globalCoordinates = problem.coords * N';

%% TODO: Generate properly distribute (i.e. not equidistant) CPs from the knot vector !!!
% In the meanwhile an easy hack is:
globalCoordinates = parametricCoordinate * 0.001;

end

