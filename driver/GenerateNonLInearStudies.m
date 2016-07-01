%% Generate Convergence Study
% Read data from .txt file and plot the convergnce study of h-FEM and
% XPODFEM w.r.t. Reference overkilled solution.

figure(10000000)
% Create axes
axes1 = axes;

referenceFile = fopen('myNonLinearReferenceResultsFile.txt','r');
formatSpec = '%f;\n';
% referenceTemperatureSolution = fscanf(referenceFile, formatSpec);
referenceTemperatureSolution = dlmread('myNonLinearReferenceResultsFile.txt');

%% hFEM plot
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for classical FEM discretization

dofsVector=[41, 81, 161, 321, 641, 1281, 2561, 5121]; % refDepth = 1,..,8
relError = [];
for depth=1:6
    formatSpec = 'myFEMNonLinearResultsFile_%d.txt';
    FEMFileID = sprintf(formatSpec,depth);
    FEMFile = fopen(FEMFileID, 'r');
    FEMTemperatureSolution = dlmread(FEMFileID);
    
    err = 0.0;
    for i=1:size(FEMTemperatureSolution,1)
        for j=1:size(FEMTemperatureSolution,2)-1
            err = err + sqrt(( referenceTemperatureSolution(i,j) - FEMTemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j)^2 );
        end
    end
    
    relError = [relError err/size(FEMTemperatureSolution,1)/(size(FEMTemperatureSolution,2)-1)];
    
end
% Create multiple lines using matrix input to plot
loglog(dofsVector(1:6),relError,'-kd');
hold on

%% X-PODFEM plot
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for XPODFEM discretization

dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
relErrorPOD = [];
pointwiseRelError = [];


for modes=1:5
    formatSpec = 'myPODXFEMNonLinearResultsFile_%d.txt';
    XFEMFileID = sprintf(formatSpec,modes);
    XFEMFile = fopen(XFEMFileID, 'r');
    XFEMTemperatureSolution = dlmread(XFEMFileID);
    
    err = 0.0;
    for i=1:size(XFEMTemperatureSolution,1)
        
        timeStepErrorAtPoint = 0.0;
        for j=2:size(XFEMTemperatureSolution,2)-1
            pointError = sqrt(( referenceTemperatureSolution(i,j-1) - XFEMTemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j-1)^2 );
            timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
        end
        pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XFEMTemperatureSolution,2)-2)];
        err = err + timeStepErrorAtPoint;
    end
    
    relErrorPOD = [relErrorPOD err/size(XFEMTemperatureSolution, 1)/(size(XFEMTemperatureSolution,2)-2)];
    
end

% Create multiple lines using matrix input to plot
loglog(dofsVectorPOD(1:5),relErrorPOD, '--bd');
hold on

% Title
title({'\centering{\quad \quad \quad Non-linear X-PODFEM}'; '\centering{using refined hFEM (depth 8) as reference}'}...
    , 'Interpreter','latex');

% Legend
legend('FEM', 'POD', 'Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on average temperature (last layer)');
xlabel('\fontname{Latin Modern Math} DOFs');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');


hold off

figure(123456789)
plot(linspace(0.0, 0.006, numel(pointwiseRelError)), pointwiseRelError);
% Title
title({'Pointwise relative error'}...
    , 'Interpreter','latex');

% % Legend
% legend('FEM', 'POD', 'Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on average temperature (last layer)');
xlabel('\fontname{Latin Modern Math} bar length');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');



figure(99999)
formatSpec = 'myPODXFEMNonLinearTimeFile_1.txt';
XFEMTimeFileID = sprintf(formatSpec,modes);
XFEMTime = dlmread(XFEMTimeFileID);
t = linspace(0, numel(XFEMTime), numel(XFEMTime)-1);
plot(t, XFEMTime(1:end-1),'-r*')
hold on

formatSpec = 'myFEMNonLinearTimeFile_4.txt';
FEMTimeFileID = sprintf(formatSpec,modes);
FEMTime = dlmread(FEMTimeFileID);
t = linspace(0, numel(FEMTime), numel(FEMTime)-1);
plot(t, FEMTime(1:end-1),'-kd')
hold on

% Title
title({'CPU Time at each time step'}...
    , 'Interpreter','latex');

% Legend
legend('POD', 'FEM', 'Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} CPU Time [msec]');
xlabel('\fontname{Latin Modern Math} Time Step');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
hold off

