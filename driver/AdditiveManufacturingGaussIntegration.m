%% 1D Multi-layer Linear Additive Manufacturing with enriched local functions
% evaluate a 1D poisson transient linear problem using X-FEM

clear all;
clc;

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
k = 56.0;                                                   % thermal conductivity [W/(m°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
heatCapacity= rho*c;                                        % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                           % source temperature [°C]


tEnd = 0.1;                                                % total time [sec]
xEnd = 0.001;                                               % length of the bar 0.001 [m]

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
nuemannRightBC = 0.0e+08;
bodySource = 0.0e+13;

for modes = 1:5
    % maxTrainingTimeSteps = timeSteps*0.5;
    % relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
    % modes = zeros(maxTrainingTimeSteps-3,1);
    % tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);
    
    myVideoSpec = 'AdditiveManufacturingPODBasis_%d.avi';
    myVideo = sprintf(myVideoSpec,modes);
    writerObj = VideoWriter(myVideo);
    writerObj.Quality = 100;
    writerObj.FrameRate = 5;
    open(writerObj);
    
    numberOfLayers = 20;
    trainingTimeSteps = 10;
    numberOfTimeStepsPerLayer = 5;

    % for trainingTimeSteps = 3:maxTrainingTimeSteps
    
    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    refinementDepth = 2;
    numberOfPODModes = modes;
    
    integrationOrder = (modes-1)^refinementDepth + 1;
    
    bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
    numberOfTimeStepsPerLayer, refinementDepth, bodySource);

    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);                        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                                 % layer spatial discretization X
    
    x_PostProcess = linspace(0.0, xEnd, (numberOfElementsInX +1) * 2.^8);                       % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, CPUTime] = backwardEulerGaussIntegrationSolver(x, x_PostProcess, bodyLoad, T0,...
        dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, refinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
        numberOfLayers, numberOfPODModes, integrationOrder);
    
    %% Post-Process
    figure(modes)
    
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
        xlim(axes1,[0 0.001]);
        ylim(axes1,[-100 3000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    hold off
    
    % Write results to a file
    formatSpec = 'myPODXFEMResultsFile_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-j));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'myPODXFEMTimeFile_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', CPUTime(i));
    end
    
    fclose(resultFile);
    
%     figure(modes + 1000)
%     % Create axes
%     axes1 = axes;
%     
%     F(size(t,2)) = struct('cdata',[],'colormap',[]);
%     
%     for i=1:size(t,2)
%         plot(x_PostProcess,temperatureSolution(:,i)')
%         
%         % Title
%         myTitleSpec = 'Temperature evolution of the bar using %d. modes';
%         myTitle = sprintf(myTitleSpec,modes);
%         title( myTitle, 'Interpreter','latex');
%         
%         % Create ylabel
%         ylabel('\fontname{Latin Modern Math} Temperature [°C]');
%         xlabel('\fontname{Latin Modern Math} length [m]');
%         
%         % Axis limit
%         xlim(axes1,[0 0.001]);
%         ylim(axes1,[-100 3000]);
%         
%         box(axes1,'on');
%         % Set the remaining axes properties
%         set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
%         
%         grid on
%         drawnow
%         F(i) = getframe;
%         writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
%     end
%     
%     close(writerObj);
    
    
    
end % modes loop


