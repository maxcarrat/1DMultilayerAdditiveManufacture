%% Poisson Transient Problem Test 1D
% evaluate a 1D poisson transient problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('XFEMMultiElement_test.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 1.0 ;                                         % density [kg/m^3]
c = 1.0;                                            % specific heat [J/(kg°C)]
k = 1.0;                                            % thermal conductivity [W/(m°C)]
T0 = 0.0;                                           % initial temperature [°C]
heatCapacity= rho*c;                                % heat capacity [kJ / kg °C]

tEnd = 8;
xEnd = 0.1e-00;

dirichletLeftBC = @(t) 0.0;
dirichletRightBC = @(t) 0.0;
rhs = 0.0;

tolerance = 1.0e-05;
maxIteration = 100;

initialNumberOfModes = 1;
maxNumberOfModes = 1;
DOFs = zeros(maxNumberOfModes, 1);

%%
% 
%  Using 1 or 2 POD modes the FEM solution and the X-PODFEM soluton
%  overlap each others. If the number of POD modes increases, e.g. using 3
%  POD modes the solution starts to oscillates. In particular we can
%  observe this phenomena around the right end of the bar. With a number of
%  POD modes >= 4 the X-PODFEM solution explodes, i.e. the matrix of the
%  eXtended system is "singular to working precision". In case of 3 POD
%  modes the size of oscillations is strongly affected by the penalty term
%  used to constrain the system, in fact the RCOND of the system matrix
%  drops from 10^-34 to 10^-24 symply decreasing the penalty term from
%  10^20 to 10^12.
%

for modes = initialNumberOfModes:maxNumberOfModes
    
    numberOfLayers = 1;
    trainingTimeSteps = 4;
    numberOfTimeStepsPerLayer = 16;                     % total time per layer 0.225 sec
    numberOfHeatingTimeSteps = 4;                       % heating laser time per layer 0.090 sec
    refinementDepth = 5;

    PODRefinementDepth = 2;

    integrationOrder = 2;
    integrationModesOrder = modes + 1; 
    
    numberOfRefinedElementsToBeKept = 1;

    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    numberOfPODModes = modes;
    
    bodyLoad = @(x, t) T0 + 1.0e+06 * x * abs(sin( pi * t));
    
    k = @(x, t, T)  k;                                                                          % thermal conductivity [W / (m K)]
    heatCapacity= @(x, T, T_last)  heatCapacity;                                                % heat capacity [kJ / kg °C]

    t = linspace(0, tEnd, numberOfTimeStepsPerLayer * numberOfLayers);                          % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                                 % layer spatial discretization X
    
    x_PostProcess = linspace(0.0, xEnd, ( numberOfElementsInX + 1 ) * 2.^5);                    % post-processing coordinates
    [X, T] = meshgrid(x_PostProcess, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, CPUTime, DOFs(modes)] = multiPhaseBackwardEulerSolverSingleLayer...
        ( x, x_PostProcess, bodyLoad, T0,...
        dirichletLeftBC, dirichletRightBC, rhs, k, heatCapacity, t, tolerance, maxIteration, numberOfRefinedElementsToBeKept,...
        refinementDepth, PODRefinementDepth, trainingTimeSteps, numberOfPODModes, integrationOrder, integrationModesOrder );
    
    %% Post-Process
    figure(modes + 111111)
    
    % Create axes
    axes1 = axes;
    
    for i=1:trainingTimeSteps
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,temperatureSolution(:,i)', 'k');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[-100 200]);
        
        box(axes1,'on');
        
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        grid on
        hold on
        
    end
    
    for i=trainingTimeSteps+1:size(temperatureSolution, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,temperatureSolution(:,i)', 'r');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[-100 200]);
        
        box(axes1,'on');
        
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        grid on
        hold on
        
    end
    
    hold off
    
   figure(modes+9999)
    % Create axes
    axes2 = axes;
    
    for i=1:trainingTimeSteps
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess, heatFlux(:,i)', 'k');
        
        % Title
        title('Heat flux distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-1.0e+04 1.0e+04]);

        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    for i=trainingTimeSteps+1:size(heatFlux, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess, heatFlux(:,i)', 'r');
        
        % Title
        title('Heat flux distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-1.0e+04 1.0e+04]);

        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    % Create ylabel
    
    
    hold off
    
    % Write results to a file
    formatSpec = 'myXFEMMultiElement_test_temperatures.txt';
    filename = sprintf(formatSpec, modes, PODRefinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write fluxes to a file
    formatSpec = 'myXFEMMultiElement_test_fluxes.txt';
    filename = sprintf(formatSpec, modes, PODRefinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(heatFlux, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', heatFlux(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'myXFEMMultiElement_test_CPUTime.txt';
    filename = sprintf(formatSpec, modes, PODRefinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', CPUTime(i));
    end
    
    fclose(resultFile);
    
    
    % Write DOFs to a file
    formatSpec = 'myXFEMMultiElement_test_DOFs.txt';
    filename = sprintf(formatSpec, modes, PODRefinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(DOFs)
        fprintf(resultFile, '%d, %\t', DOFs(i));
    end
    
    fclose(resultFile);

end % modes loop