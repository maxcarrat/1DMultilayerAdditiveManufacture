function heatFlux = evaluateHeatFlux( problem, temperature)
% heatFlux = evaluateHeatFlux(problem, temperature) evaluates the numerical solution
%   problem = struct that defines the boundary value problem
%   temeperature = solutions of the problem    

    heatFlux = zeros(size(problem.coords(1), 1));

    for i = 1:size(problem.coords,1)
        heatFlux(i) = problem.k * problem.basis_fun(problem.coords(i), 1, 1) * temperature(i);
    end

end
