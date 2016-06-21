%% 1D Multi-layer Linear Additive Manufacturing 
% evaluate a 1D poisson transient linear problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('AdditiveManufacturingHFEM.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 1000.0;                                          % density [kg/m^3]
c = 1.0;                                               % specific heat [J/(kg°C)]
k = 1.0;                                               % thermal conductivity [W/(m°C)]
T0 = 20.0;                                             % Initial temperature [°C]
heatCapacity= rho*c;                                   % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                      % source temperature [°C]

tEnd = 10.0e-03;
xEnd = 0.1;

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
rhs = @(x, t) 0.0;%1.0e+08;

for depth = 6:6
    % maxTrainingTimeSteps = timeSteps*0.5;
    % relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
    % modes = zeros(maxTrainingTimeSteps-3,1);
    % tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);
    
    numberOfLayers = 20;
    trainingTimeSteps = 20;
    numberOfTimeStepsPerLayer = 10;
    
    integrationOrder = 2;
    
    % for trainingTimeSteps = 3:maxTrainingTimeSteps
    
    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    refinementDepth = depth;
    numberOfPODModes = 0;
    
    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                 % layer spatial discretization X
    
    x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref
    x_PostProcess = linspace(0.0, xEnd, (numberOfElementsInX +1) * 2^8);          % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    [X_ref, T_ref] = meshgrid(x_ref, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, modes(trainingTimeSteps-2)] = backwardEulerGaussIntegrationSolver(x, x_PostProcess, rhs, T0,...
        dirichletLeftBC, dirichletRightBC, k, heatCapacity, t, refinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
        numberOfLayers, numberOfPODModes, integrationOrder);
    
    %% Post-Process
    overkilledInternalEnergy = 4.221615373269894e+07;
    
    relErrorEnergyNorm(trainingTimeSteps-2) = sqrt((internalEnergy(end)-overkilledInternalEnergy)^2)/sqrt(overkilledInternalEnergy^2);
    
    figure(depth+6)
    
    % Create axes
    axes1 = axes;
    
    for i=1:size(temperatureSolution, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,temperatureSolution(:,i)');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 0.1]);
        ylim(axes1,[-1000 5000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    % Create ylabel
    
    
    hold off
    
    % Write results to a file
    formatSpec = 'myFEMResultsFile_%d.txt';
    filename = sprintf(formatSpec,depth);
    
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
%         for j=1:size(temperatureSolution, 2)
            fprintf(resultFile, '%d;\n', temperatureSolution(i, end-1));
%         end
%         fprintf(resultFile, ';\n');
    end
    fclose(resultFile);
    
    
    
end % depth loop


figure(10101)

% Create axes
axes1 = axes;

F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_PostProcess,temperatureSolution(:,i)')
    
    % Title
    title('Temperature evolution of the bar', 'Interpreter','latex');
    
    % Create ylabel
    ylabel('\fontname{Latin Modern Math} Temperature [°C]');
    xlabel('\fontname{Latin Modern Math} length [m]');
    
    % Axis limit
    xlim(axes1,[0 0.1]);
    ylim(axes1,[-1000 5000]);
    
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
    
    grid on
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);