function heatFlux = evaluateHeatFlux( problem, temperature, timeStep)
% heatFlux = evaluateHeatFlux(problem, temperature) evaluates the numerical solution
%   problem = struct that defines the boundary value problem
%   temeperature = solutions of the problem    

    coords = problem.coords;
    heatFlux = zeros(size(coords, 2), 1);
   
    gradientOperator =@(x) getGradientOperator(problem, x);
    B_0 = gradientOperator(mapGlobalToLocal(coords(1), -1, 1));
    
    heatFlux(1) = problem.k * B_0(1) * temperature(1, timeStep)';

    for i=2:size(coords,2)
        B_i = gradientOperator(mapGlobalToLocal(coords(i)-coords(i-1), -1, 1));
        heatFlux(i) = problem.k * B_i(i-1) * temperature(i, timeStep)';
    end
    
end
