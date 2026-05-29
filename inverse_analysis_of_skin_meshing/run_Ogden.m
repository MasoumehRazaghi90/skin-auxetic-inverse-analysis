function results = run_Ogden()

%%
% Plot settings
fontSize=20;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=15;
markerSize2=40;
markerSize3=50;
lineWidth1=3;
lineWidth2=3;
cMap=viridis(20);
numPlotPointsGraphs = 100;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z
febioLogFileName_stretch=[febioFebFileNamePart,'_stretch_out.txt']; %Log file name for exporting stretch U_z

% Define data paths
loadPath_experimental = fullfile(defaultFolder,'data','Ni_Annaidh_2012');
dataName_1 = fullfile(loadPath_experimental,'Ni_Annaidh_stress_stretch_perp.mat');
dataName_2 = fullfile(loadPath_experimental,'Ni_Annaidh_stress_stretch_para.mat');

%Specifying dimensions and number of elements
sampleWidth=6;
sampleThickness=2.25;
sampleHeight=33;
% pointSpacings=10*ones(1,3);
pointSpacings = [1 0.5 3];
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));

%True material parameter set
k_factor=100;
c1_true=0.5;
m1_true=6;
k_true=c1_true*k_factor;

%Initial material parameter set
c1_ini=c1_true;
m1_ini=m1_true;
c2_ini = 0.05;    % initial guess for second term
m2_ini = 12;   % initial guess (negative m)
k_ini=c1_ini*k_factor;

evalMode = 2; % 1= test, 2=optimise

P = [c1_ini, m1_ini, c2_ini, m2_ini];


% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%Optimisation settings
maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';
optimisationMethod = 2; % 1= fminsearch, Nelder-Mead, 2=lsqnonlin, Levenberg-Marquart


%% LOAD EXPERIMENTAL DATA

data_perp = load(dataName_1);
data_para = load(dataName_2);

stretch_exp_perp_raw = data_perp.X1;
stress_exp_perp_raw = data_perp.Y1;

stretch_exp_para_raw = data_para.X1;
stress_exp_para_raw = data_para.Y1;

%% --- Make ONE mean experimental curve up to max stretch of perp ---

stretchLoad = min(max(stretch_exp_perp_raw), max(stretch_exp_para_raw));
stretch_common = linspace(1, stretchLoad, numPlotPointsGraphs);

stress_perp_on_common = interp1(stretch_exp_perp_raw, stress_exp_perp_raw, ...
    stretch_common, 'pchip');

stress_para_on_common = interp1(stretch_exp_para_raw, stress_exp_para_raw, ...
    stretch_common, 'pchip');

% keep parallel stress constant after its last point
lastStress = stress_exp_para_raw(end);
stress_para_on_common(stretch_common > max(stretch_exp_para_raw)) = lastStress;

stress_mean = 0.5*(stress_perp_on_common + stress_para_on_common);

stretch_exp_mean_raw = stretch_common(:);
stress_exp_mean_raw  = stress_mean(:);

% Use ONLY perp max displacement (since stretch limit is perp)
stretchLoad_mean = stretchLoad;
displacementMagnitude_mean = (stretchLoad_mean*sampleHeight) - sampleHeight;

%%
%% Interpolated mean experimental curve + FEA

stretch_plot_mean = linspace(1, stretchLoad_mean, numPlotPointsGraphs);

stress_plot_mean = interp1(stretch_exp_mean_raw, stress_exp_mean_raw,...
    stretch_plot_mean,'pchip');

hf1 = cFigure; hold on;
title('Stretch stress curve optimized','FontSize',fontSize);

xlabel('\lambda Stretch [.]','FontSize',fontSize);
ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize);
zlabel('Z','FontSize',fontSize); 
hold on;

% Experimental mean data
Hn(1) = plot(stretch_exp_mean_raw,stress_exp_mean_raw,'k.','MarkerSize',markerSize2);

% Interpolated experimental mean curve
Hn(2) = plot(stretch_plot_mean,stress_plot_mean,'k-','LineWidth',lineWidth1);

% FEA curve (updated during optimisation)
Hn(3) = plot(NaN,NaN,'bx','MarkerSize',markerSize2,'LineWidth',1.5);

legend(Hn,{'Ni Annaidh Experimental mean','Interpolated mean','FEA mean'},'Location','northwest');

view(2);
axis tight;
axis square;
box on;
grid on;
set(gca,'FontSize',fontSize);
drawnow;

%% PLOT: show only MEAN experimental curve (one target curve)

stretch_plot = stretch_exp_mean_raw(:);
stress_plot  = stress_exp_mean_raw(:);

% hf1 = cFigure; hold on;
% title('Stretch-stress (mean curve) + FEA','FontSize',fontSize);
% xlabel('\lambda Stretch [.]','FontSize',fontSize);
% ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize);
% 
% % Experimental mean curve
% Hn(1) = plot(stretch_plot, stress_plot, 'k-', 'LineWidth', lineWidth1);
% 
% % FEA curve (will be updated inside objective function)
% Hn(2) = plot(NaN, NaN, 'b-', 'LineWidth', lineWidth1);
% 
% legend(Hn, {'Experimental mean', 'FEA'}, 'Location','northwest');
% view(2); axis tight; axis square; box on; grid on;
% set(gca,'FontSize',fontSize);
% drawnow;
%% CREATING MESHED BOX

%Create box 1
boxDim=[sampleWidth sampleThickness sampleHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
[box1]=hexMeshBox(boxDim,boxEl);
E=box1.E;
V=box1.V;
Fb=box1.Fb;
faceBoundaryMarker=box1.faceBoundaryMarker;

X=V(:,1); Y=V(:,2); Z=V(:,3);
VE=[mean(X(E),2) mean(Y(E),2) mean(Z(E),2)];

elementMaterialIndices=ones(size(E,1),1);

%%

% Plotting boundary surfaces

cFigure; hold on;
title('Model surfaces','FontSize',fontSize);
gpatch(Fb,V,faceBoundaryMarker,'k',0.5);
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);
drawnow;

%% DEFINE BC's

%Define supported node sets
logicFace=faceBoundaryMarker==1;
Fr=Fb(logicFace,:);
bcSupportList_X=unique(Fr(:));

logicFace=faceBoundaryMarker==3;
Fr=Fb(logicFace,:);
bcSupportList_Y=unique(Fr(:));

logicFace=faceBoundaryMarker==5;
Fr=Fb(logicFace,:);
bcSupportList_Z=unique(Fr(:));

%Prescribed displacement nodes
logicPrescribe=faceBoundaryMarker==6;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));

%%
% Visualize BC's
cFigure; hold on;
title('Complete model','FontSize',fontSize);

gpatch(Fb,V,'kw','k',0.5);
plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize1);
plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize1);

axisGeom(gca,fontSize);
drawnow;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings
[febio_spec]=febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version='4.0';

%Module section
febio_spec.Module.ATTR.type='solid';

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1_ini;
febio_spec.Material.material{1}.m1=m1_ini;
febio_spec.Material.material{1}.c2=c2_ini;
febio_spec.Material.material{1}.m2=m2_ini;
febio_spec.Material.material{1}.k=k_ini;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList_X';
nodeSetName2='bcSupportList_Y';
nodeSetName3='bcSupportList_Z';
nodeSetName4='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList_X);

febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList_Y);

febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcSupportList_Z);

febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
febio_spec.Mesh.NodeSet{4}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_x';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=0;
febio_spec.Boundary.bc{1}.z_dof=0;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_y';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof=0;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_z';
febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
febio_spec.Boundary.bc{3}.x_dof=0;
febio_spec.Boundary.bc{3}.y_dof=0;
febio_spec.Boundary.bc{3}.z_dof=1;

febio_spec.Boundary.bc{4}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{4}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
febio_spec.Boundary.bc{4}.dof='z';
febio_spec.Boundary.bc{4}.value.ATTR.lc=1;
% febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitude_perp;
febio_spec.Boundary.bc{4}.value.VAL = displacementMagnitude_mean;
febio_spec.Boundary.bc{4}.relative=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stretch;
febio_spec.Output.logfile.element_data{2}.ATTR.data='Uz';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';

% Plotfile section
febio_spec.Output.plotfile.compression=0;

%% Creating febio analysis structure

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=0; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

%% Create structures for optimization

%What should be known to the objective function:
objectiveStruct.stretch_exp_mean_raw = stretch_exp_mean_raw;
objectiveStruct.stress_exp_mean_raw  = stress_exp_mean_raw;

objectiveStruct.displacementMagnitude_mean = displacementMagnitude_mean;

objectiveStruct.febioAnalysis=febioAnalysis;
objectiveStruct.febio_spec=febio_spec;
objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.febioLogFileName_disp = fullfile(savePath,febioLogFileName_disp);
objectiveStruct.febioLogFileName_stress = fullfile(savePath,febioLogFileName_stress);
objectiveStruct.febioLogFileName_stretch = fullfile(savePath,febioLogFileName_stretch);


objectiveStruct.parNormFactors=P; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=P; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim = [
    1e-4   10;   % c1
    1      100;     % m1
    1e-4   10;   % c2
    1      100    % m2
    ];

objectiveStruct.k_factor=k_factor;


objectiveStruct.hf = hf1;

objectiveStruct.h  = Hn(3);

objectiveStruct.method = optimisationMethod;

Pn=P./objectiveStruct.parNormFactors;
P_opt = P;   % default if no optimisation


switch evalMode
    case 1
        [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct);
        SSE = sum(errorVal.^2);
        n = numel(errorVal);
        stdError = sqrt(SSE/(n-1));
        % stretch_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
        % stress_sim  = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));


    case 2
        %% start optimization
        switch optimisationMethod
            case 1 %fminsearch and Nelder-Mead
                OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
                OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
            case 2 %lsqnonlin and Levenberg-Marquardt
                OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);
        end

        %% Unnormalize and constrain parameters

        [errorVal,simData] = objectiveFunctionIFEA(Pn_opt, objectiveStruct);
        SSE = sum(errorVal.^2);
        n = numel(errorVal);
        stdError = sqrt(SSE/(n-1));

        % stretch_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
        % stress_sim  = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));

        P_opt=Pn_opt.*objectiveStruct.parNormFactors; %Scale back, undo normalization

        %Constraining parameters
        for q=1:1:numel(P_opt)
            [P_opt(q)]=boxconstrain(P_opt(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
        end

        disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
        disp(['P_opt=',disp_text]);

end

%% Import FEBio results

% Displacements
N_disp_mat1=simData(1).dataStruct_disp.data; %Displacement
timeVec1=simData(1).dataStruct_disp.time; %Time
V_DEF1=N_disp_mat1+repmat(V,[1 1 size(N_disp_mat1,3)]); %Deformed coordinate set
DN_magnitude1=sqrt(sum(N_disp_mat1(:,:,end).^2,2)); %Current displacement magnitude


%%
% % Plotting the simulated results using |anim8| to visualize and animate
% % deformations
% Animate deformation (single simulation)

hf = cFigure;
title('Displacement magnitude [mm]','Interpreter','latex')

hp = gpatch(Fb,V_DEF1(:,:,end),DN_magnitude1,'k',1,2);
hp.Marker='.';
hp.MarkerSize=markerSize2;
hp.FaceColor='interp';

gpatch(Fb,V,0.5*ones(1,3),'none',0.25);

axisGeom(gca,fontSize);
colormap(cMap); 
colorbar;

caxis([0 max(DN_magnitude1)]);
axis(axisLim(V_DEF1));

view(140,30);
camlight headlight;

animStruct = struct;
animStruct.Time = timeVec1;

for qt = 1:size(N_disp_mat1,3)

    DN_magnitude = sqrt(sum(N_disp_mat1(:,:,qt).^2,2));

    animStruct.Handles{qt} = hp;
    animStruct.Props{qt}   = {'Vertices','CData'};
    animStruct.Set{qt}     = {V_DEF1(:,:,qt),DN_magnitude};

end

anim8(hf,animStruct);
drawnow;



%% Stress plot (ONE simulation only)

%Access data
E_stress_mat1=simData(1).dataStruct_stress.data;
E_stretch_mat1=simData(1).dataStruct_stretch.data;
[CV1]=faceToVertexMeasure(E,V,E_stress_mat1(:,:,end));

hf = cFigure;
title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')

hpS = gpatch(Fb, V_DEF1(:,:,end), CV1, 'k', 1, 2);
hpS.Marker='.'; hpS.MarkerSize=markerSize2; hpS.FaceColor='interp';
gpatch(Fb, V, 0.5*ones(1,3), 'none', 0.25);

axisGeom(gca,fontSize);
colormap(gca,cMap); colorbar;
caxis([min(E_stress_mat1(:)) max(E_stress_mat1(:))]);
axis(axisLim(V_DEF1));
view(140,30);
camlight headlight;

animStruct = struct;
animStruct.Time = timeVec1;

for qt = 1:size(N_disp_mat1,3)
    CV1t = faceToVertexMeasure(E, V, E_stress_mat1(:,:,qt));
    animStruct.Handles{qt} = hpS;
    animStruct.Props{qt}   = {'Vertices','CData'};
    animStruct.Set{qt}     = {V_DEF1(:,:,qt), CV1t};
end

anim8(hf, animStruct);
drawnow;
%%

function [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct)

%% Access input structure

% --- read mean experimental data ---
stretch_exp = objectiveStruct.stretch_exp_mean_raw(:);
stress_exp  = objectiveStruct.stress_exp_mean_raw(:);

% --- read FE model info ---
displacementMagnitude_mean = objectiveStruct.displacementMagnitude_mean;

febioAnalysis = objectiveStruct.febioAnalysis;
febio_spec = objectiveStruct.febio_spec;
febioFebFileName = objectiveStruct.febioFebFileName;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
febioLogFileName_stress = objectiveStruct.febioLogFileName_stress;
febioLogFileName_stretch = objectiveStruct.febioLogFileName_stretch;

parNormFactors = objectiveStruct.parNormFactors;
Pb_struct = objectiveStruct.Pb_struct;
k_factor = objectiveStruct.k_factor;

hf = objectiveStruct.hf;
h = objectiveStruct.h;

%% Set/get material parameters
% Unnormalize and constrain parameters
P = Pn.*parNormFactors; %Scale back, undo normalization

%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),Pb_struct.xxlim(q,1),Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
end

disp(['Trying   : ',sprintf(repmat('%6.8e ',[1,numel(P)]),P)]);

c1 = P(1);
m1 = P(2);
c2 = P(3);
m2 = P(4);

% 2-term Ogden

% Set material parameters
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c2;
febio_spec.Material.material{1}.m2=m2;
febio_spec.Material.material{1}.k=c1 * k_factor;



nElem = size(febio_spec.Mesh.Elements{1}.elem.VAL,1);


febio_spec.Boundary.bc{4}.value.VAL = displacementMagnitude_mean;

febioStruct2xml(febio_spec, febioFebFileName);
runFlag = runMonitorFEBio(febioAnalysis);

if ~runFlag
    simData(1).dataStruct_disp    = NaN;
    simData(1).dataStruct_stress  = NaN;
    simData(1).dataStruct_stretch = NaN;
    errorVal = NaN(size(stress_exp));
    return
end


simData(1).dataStruct_disp    = importFEBio_logfile(febioLogFileName_disp, 0, 1);
simData(1).dataStruct_stress  = importFEBio_logfile(febioLogFileName_stress, 0, 1);
simData(1).dataStruct_stretch = importFEBio_logfile(febioLogFileName_stretch, 0, 1);



stretch_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
stress_sim = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));

stress_sim_on_exp = interp1(stretch_sim, stress_sim, stretch_exp,'pchip','extrap');

errorVal = stress_exp(2:end) - stress_sim_on_exp(2:end);

figure(hf);
set(h,'XData',stretch_sim,'YData',stress_sim);
drawnow;

n = numel(errorVal); % number of residuals


if objectiveStruct.method == 1
    SSE = sum(errorVal.^2);
    errorVal = SSE;
end

end
stretch_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
stress_sim  = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));

results.stretch_perp = stretch_sim;
results.stress_perp  = stress_sim;

results.stretch_para = stretch_sim;
results.stress_para  = stress_sim;

results.stretch_exp_perp = stretch_exp_perp_raw;
results.stress_exp_perp  = stress_exp_perp_raw;

results.stretch_exp_para = stretch_exp_para_raw;
results.stress_exp_para  = stress_exp_para_raw;

results.stretch_exp_mean = stretch_exp_mean_raw;
results.stress_exp_mean  = stress_exp_mean_raw;

results.SSE = SSE;
results.stdError = stdError;
results.P_opt = P_opt;
end
