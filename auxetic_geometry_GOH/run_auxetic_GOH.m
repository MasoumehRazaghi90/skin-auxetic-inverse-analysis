
clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8; %transparency
markerSize=40; %For plotted points
markerSize2=10; %For nodes on patches
lineWidth1=1; %For meshes
cMap=spectral(250); %colormap


%% Control parameters

meshSearchTolerance=1e-3; %Tolerance for logic on mesh searching thresholds


%Material parameter set
materialType='HGO';

c_ini     = 1.2098951004; % Initial slope
k1_ini    = 61.7903625549; % "Modulus of fibers"
k2_ini    = 3.4106088361e-7; %  "non-linearity" or J-shape
kappa_ini = 0.2985957214; % [0,1/3], "distance between graphs"
k_ini     = 1.2098951004*1000; % Bulk modulus
gamma     = 41;
LangerAngle_perp = 0/180*pi; % Angle of Langer line with loading direction
LangerAngle_para = 90/180*pi; % Angle of Langer line with loading direction
angleSelection=1;%1 perpendicular to LL and 2 parallel
LangerAngles = [LangerAngle_perp, LangerAngle_para];

% E_youngs1=1; %Material Young's modulus
% nu1=0.4; %Material Poisson's ratio


% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=20; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=2; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
runMode='internal';

% Path names
defaultFolder = fileparts(mfilename('fullpath'));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting principal stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%% geometries

% pointSpacing =2;
% distRefine = [3];

% pointSpacing =2;
% distRefine = [4 3];

% pointSpacing =1;
% distRefine = [2];

% pointSpacing =1;
% distRefine = [3 2];

% pointSpacing =2;
% distRefine = [5 4 3];

% pointSpacing =1;
% distRefine = [4 3 2];

% pointSpacing =2;
% distRefine = [6 5 4 3];

% pointSpacing =1;
% distRefine = [5 4 3 2];

% pointSpacing =2;
% distRefine = [3.5 2.5];

% pointSpacing =1;
% distRefine = [2.5 1.5];

% pointSpacing =0.5;
% distRefine = [2 1.5 1];

% pointSpacing =0.25;
% distRefine = [1.5 1 0.5];

% pointSpacing =0.5;
% distRefine = [2 1];

% pointSpacing =3;
% distRefine = [5 4];

% pointSpacing =0.5;
% distRefine = [3 2];

% pointSpacing =0.5;
% distRefine = [4 3];

pointSpacing =0.5;
distRefine = [3 1];
% 

geometryID = 5;   % choose your geometry (1–5)
geo = allGeometries(geometryID, pointSpacing, distRefine);
savePath = fullfile(defaultFolder,'data',sprintf('geometry_%d',geometryID));
if ~exist(savePath,'dir'); mkdir(savePath); end

F = geo.F;
V = geo.V;
indRefine = geo.indRefine;
V_refine = geo.V_refine;
V1 = geo.V1;
p = geo.p;
d = geo.d;
t = geo.t;

%FEBio Parameters
h0 = max(V(:,2)) - min(V(:,2));  % Initial mesh height
displacementMagnitude = 61.59/4;  % Slightly boost to avoid lambda_y = 1

%% Creating triangle Mesh

for i = 1:1:length(distRefine)
    % VF = patchCentre(F,V);
    % [distCorners,indMin] = minDist(VF,V_refine);

    %Compute distances on mesh description
    distCorners_V=meshDistMarch(F,V,indRefine);
    distCorners=vertexToFaceMeasure(F,distCorners_V);

    cFigure;
    hold on;
    % % title('Colormapped face colors');
    Hp=gpatch(F,V, distCorners_V,'k');
    Hp.FaceColor='interp';
    plotV(V(indRefine,:),'r.','Markersize',25);
    plotV(V_refine,'b.','Markersize',25);
    axisGeom;
    drawnow;


    logicRefine = distCorners<=distRefine(i);
    logicRefine = triSurfLogicSharpFix(F, logicRefine, 3);

    % [F,V,C_type,indIni,C]=subTriDual(F,V,logicRefine);
    inputStruct.F=F;
    inputStruct.V=V;
    inputStruct.indFaces=find(logicRefine);
    % inputStruct.f=f;
    [outputStruct]=subTriLocal(inputStruct);
    F=outputStruct.F;
    V=outputStruct.V;

    %Smooth with constraints
    Eb = patchBoundary(F);
    indRigid = unique(Eb(:));
    smoothParameters.n=25; %Number of iterations
    smoothParameters.Method='LAP'; %Smoothing method
    smoothParameters.RigidConstraints=indRigid; %Vertices to not keep constant
    V=patchSmooth(F,V,[],smoothParameters);
end
V = V(:,[1,2]);

cFigure;
hold on;
% title('Colormapped face colors');
gpatch(F,V,'gw','k');
gpatch(F,V,'w','k');
pointAnnotate(V1,1:1:length(V1),'fontsize',15)
axisGeom;
drawnow;


%%
F=[F;fliplr(F)+length(V)];
Vc=V;
Vc(:,1)=-Vc(:,1);
V=[V;Vc];


F=[F;fliplr(F)+length(V)];
Vc=V;
Vc(:,2)=-Vc(:,2);
V=[V;Vc];


%% replicate mesh

[FT,VT,CT]=replicatemesh(F,V,p,d);

%%

%Smooth with constraints
Eb = patchBoundary(FT);
indRigid = unique(Eb(:));
smoothParameters2.n=25; %Number of iterations
smoothParameters2.Method='LAP'; %Smoothing method
smoothParameters2.RigidConstraints=indRigid; %Vertices to not keep constant
VT=patchSmooth(FT,VT,[],smoothParameters2);

%% Visualizing triangle Mesh
cFigure;
hold on;
% title('Colormapped face colors');
gpatch(FT,VT,CT,'k');
plotV(V1,'k.-','linewidth',3,'markersize', 15)
% pointAnnotate(V1,1:1:length(V1),'fontsize',15)
axisGeom;
colormap(gjet(prod(p))); icolorbar;
drawnow;


%% Thickened the solidmesh
numSteps=ceil(t/pointSpacing);
[E,V,Fp1,Fp2]=patchThick(FT,VT,1,t,numSteps);

%Use element2patch to get patch data
F=element2patch(E,[],'penta6');



%%

%Get boundary faces (two sets due to pentahedra)
indb = tesBoundary(F); %Cell containing boundary face indices
Fb = {F{1}(indb{1},:),F{2}(indb{2},:)}; %Cell containing boundary faces

%%
% Plotting meshed model
cFigure; hold on;
title('The meshed model','FontSize',fontSize);

gpatch(F,V,'g','k',0.5); %All faces
gpatch(Fb,V,'gw','k',faceAlpha1); %Boundary faces

patchNormPlot(F,V); %Visualise normal directions

axisGeom(gca,fontSize);
camlight headlight;
drawnow;


%% Defining the boundary conditions

%Prescribed displacement nodes
bcPrescribeList=find(V(:,2)> max(V(:,2))-meshSearchTolerance);
bcSupportList=find(V(:,2)< min(V(:,2))+meshSearchTolerance);

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','k',0.5);

hl(1)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(V(bcSupportList,:),'g.','MarkerSize',markerSize);

legend(hl,{'Tensile Positive', 'Tensile Negative'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;


%% DEFINE FIBRE DIRECTIONS
nElem=length(E);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%
setenv('GIBBON_cFigure_close','0');

angles_deg = 0;  % Loop from 0 to 90 degrees

for ang = angles_deg
    % Define fiber direction by rotating vector [1 0 0] around z-axis
    R = euler2DCM([0,0,-ang*pi/180]);
    e1 = (R*[1 0 0]')';  % Rotated fiber direction

    e1_dir = repmat(e1, nElem, 1);
    e2_dir = repmat([0 0 1], nElem, 1);  % Keep this fixed
    e3_dir = cross(e1_dir, e2_dir, 2);   % Orthogonal vector

    % Plot material directions
    VE = patchCentre(E, V);  % Make sure VE is defined

    cFigure; hold on;
    title(sprintf('Fiber direction = %d°', ang), 'FontSize', 20);
    gpatch(Fb, V, 'kw', 'none', 0.25);
    quiverVec(VE, e1_dir, mean(pointSpacing), 'r');
    quiverVec(VE, e2_dir, mean(pointSpacing)/2, 'g');
    quiverVec(VE, e3_dir, mean(pointSpacing)/2, 'b');
    axisGeom(gca, fontSize);
    camlight headlight;
    drawnow;



    %%
    % Visualizing material directions

    [VE]=patchCentre(E,V);

    hf=cFigure; hold on;
    title(sprintf('Material directions  |  Fiber Angle = %d°', ang), 'FontSize', fontSize);

    gpatch(Fb,V,'kw','none',0.25);
    hf(1)=quiverVec(VE,e1_dir,mean(pointSpacing),'r');
    hf(2)=quiverVec(VE,e2_dir,mean(pointSpacing)/2,'g');
    hf(3)=quiverVec(VE,e3_dir,mean(pointSpacing)/2,'b');

    legend(hf,{'e1-direction','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;


    legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;
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
    switch materialType
        case 'HGO'

            febio_spec.Material.material{1}.ATTR.type = 'uncoupled solid mixture';
            materialName1='Material1_HGO';
            febio_spec.Material.material{1}.ATTR.name=materialName1;
            febio_spec.Material.material{1}.ATTR.type = 'solid mixture';
            febio_spec.Material.material{1}.ATTR.id=1;
            
            
            % GOH‑specific material parameters
            febio_spec.Material.material{1}.solid{1}.ATTR.type='HGO unconstrained';
            febio_spec.Material.material{1}.solid{1}.c     = c_ini/2;
            febio_spec.Material.material{1}.solid{1}.k1    = k1_ini;
            febio_spec.Material.material{1}.solid{1}.k2    = k2_ini;
            febio_spec.Material.material{1}.solid{1}.gamma = gamma;
            febio_spec.Material.material{1}.solid{1}.kappa = kappa_ini;
            febio_spec.Material.material{1}.solid{1}.k     = k_ini;

            
            % GOH‑specific material parameters
            febio_spec.Material.material{1}.solid{1}.ATTR.type='HGO unconstrained';
            febio_spec.Material.material{1}.solid{1}.c     = c_ini/2;
            febio_spec.Material.material{1}.solid{1}.k1    = k1_ini;
            febio_spec.Material.material{1}.solid{1}.k2    = k2_ini;
            febio_spec.Material.material{1}.solid{1}.gamma = -gamma;
            febio_spec.Material.material{1}.solid{1}.kappa = kappa_ini;
            febio_spec.Material.material{1}.solid{1}.k     = k_ini;


        case 'neo-Hookean'
            materialName1='Material1_neo-Hookean';
            febio_spec.Material.material{1}.ATTR.name=materialName1;
            febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
            febio_spec.Material.material{1}.ATTR.id=1;
            febio_spec.Material.material{1}.E=E_youngs1;
            febio_spec.Material.material{1}.v=nu1;
    end
    % Mesh section
    % -> Nodes

    %%Area of interest
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

    % -> Elements
    partName1='KirigamiGripper';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type='penta6'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

    % -> NodeSets
    nodeSetName1='bcPrescribeList1';
    nodeSetName2='bcSupportList2';
    nodeSetName3='bcSupportList3';

    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList);

    febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
    febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList);

    febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
    febio_spec.Mesh.NodeSet{3}.VAL=bcSupportList(1);

    %MeshDomains section
    febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
    febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

    %MeshData section
    % -> ElementData
    febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
    febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

    for q=1:1:size(E,1)
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir(q,:);
    end
    %Boundary condition section
    % -> Fix boundary conditions
    febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xz';
    febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.x_dof=0;
    febio_spec.Boundary.bc{1}.y_dof=0;
    febio_spec.Boundary.bc{1}.z_dof=1;

    febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement';
    febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
    febio_spec.Boundary.bc{2}.x_dof=0;
    febio_spec.Boundary.bc{2}.y_dof=1;
    febio_spec.Boundary.bc{2}.z_dof=1;

    febio_spec.Boundary.bc{3}.ATTR.name='prescibed_displacement_y';
    febio_spec.Boundary.bc{3}.ATTR.type='prescribed displacement';
    febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{3}.dof='y';
    febio_spec.Boundary.bc{3}.value.ATTR.lc=1;
    febio_spec.Boundary.bc{3}.value.VAL=displacementMagnitude;
    febio_spec.Boundary.bc{3}.relative=0;

    febio_spec.Boundary.bc{4}.ATTR.name='zero_displacement';
    febio_spec.Boundary.bc{4}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName3;
    febio_spec.Boundary.bc{4}.x_dof=1;
    febio_spec.Boundary.bc{4}.y_dof=0;
    febio_spec.Boundary.bc{4}.z_dof=0;


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
    logNamePrefix = sprintf('geom%d_ang_%03d',geometryID,ang);
    febioLogFileName        = [logNamePrefix, '.txt'];
    febioLogFileName_disp   = [logNamePrefix, '_disp_out.txt'];
    febioLogFileName_stress_prin = [logNamePrefix, '_stress_prin_out.txt'];
    febioLogFileName_force  = [logNamePrefix, '_force_out.txt'];

    febio_spec.Output.logfile.ATTR.file = febioLogFileName;

    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

    febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
    febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
    febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';


    % Plotfile section
    febio_spec.Output.plotfile.compression=0;

    %% Quick viewing of the FEBio input file structure
    % The |febView| function can be used to view the xml structure in a MATLAB
    % figure window.

    %%
    %%|febView(febio_spec); %Viewing the febio file|

    %% Exporting the FEBio input file
    % Exporting the febio_spec structure to an FEBio input file is done using
    % the |febioStruct2xml| function.

    febFileName = fullfile(savePath, sprintf('%s_ang_%03d.feb', febioFebFileNamePart, ang));
    febioStruct2xml(febio_spec, febFileName);
 %Exporting to file and domNode
    %system(['gedit ',febioFebFileName,' &']);

    %% Running the FEBio analysis
    % To run the analysis defined by the created FEBio input file the
    % |runMonitorFEBio| function is used. The input for this function is a
    % structure defining job settings e.g. the FEBio input file name. The
    % optional output runFlag informs the user if the analysis was run
    % succesfully.

    febioAnalysis.run_filename=febFileName; %The input file name
    febioAnalysis.run_logname=fullfile(savePath,febioLogFileName); %The name for the log file
    febioAnalysis.disp_on=0; %Display information on the command window
    febioAnalysis.runMode=runMode;
    febioAnalysis.maxLogCheckTime=10; %Max log file checking time

    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    %% Import FEBio results

    if runFlag==1 %i.e. a succesful run
        %%
        % Importing nodal forces from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),0,1);

        %Access data
        N_force_mat=dataStruct.data; %Displacement
        timeVec=dataStruct.time; %Time

        %%
        % Importing nodal displacements from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

        %Access data
        N_disp_mat=dataStruct.data; %Displacement

        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

        %%
        %% True Poisson's Ratio vs. Axial Stretch

        % Reference dimensions
        w0 = max(V(:,1)) - min(V(:,1));  % Initial width (x-direction)
        h0 = max(V(:,2)) - min(V(:,2));  % Initial height (y-direction)
        A0 = w0 * h0;                    % Initial area

        nSteps = size(V_DEF,3);
        lambda_x_vec = zeros(nSteps,1);
        lambda_y_vec = zeros(nSteps,1);
        nu_true      = zeros(nSteps,1);
        J_vec        = zeros(nSteps,1);

        for i = 1:nSteps
            V_temp = V_DEF(:,:,i);

            w = max(V_temp(:,1)) - min(V_temp(:,1)); % Current width
            h = max(V_temp(:,2)) - min(V_temp(:,2)); % Current height
            A = w * h;

            lambda_x_vec(i) = w / w0;
            lambda_y_vec(i) = h / h0;
            J_vec(i)        = A / A0;%area ratio

            % True Poisson’s ratio (log-based definition)
            if lambda_y_vec(i) > 0 && lambda_x_vec(i) > 0
                nu_true(i) = -log(lambda_x_vec(i)) / log(lambda_y_vec(i));
            else
                nu_true(i) = NaN;
            end
        end

        % Plot true Poisson’s ratio
        cFigure; hold on;
        plot(lambda_y_vec, nu_true, 'b-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('\nu_{\rm true} (Poisson''s Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Poisson''s Ratio vs. Stretch  |  Fiber Angle = %d°', ang), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;
        %% Area Ratio J vs. Axial Stretch

        cFigure; hold on;
        plot(lambda_y_vec, J_vec, 'k-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('J (Area Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Area Ratio J vs. Axial Stretch  |  Fiber Angle = %d°', ang), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;

        %%
        % Importing element stress from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),0,1);

        %Access data
        E_stress_mat=dataStruct.data;

        E_stress_mat_VM=sqrt(( (E_stress_mat(:,1,:)-E_stress_mat(:,2,:)).^2 + ...
            (E_stress_mat(:,2,:)-E_stress_mat(:,3,:)).^2 + ...
            (E_stress_mat(:,1,:)-E_stress_mat(:,3,:)).^2  )/2); %Von Mises stress

        %%
        % Compute element area (only once)
        A = patchArea(E, V_DEF(:,:,1));
        A_tot = sum(A);

        % Storage for p90 values for each time step
        p90_vec = zeros(size(N_disp_mat,3),1);

        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,end));

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy

        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('$\sigma_{vm}$ [MPa]','Interpreter','Latex')
        hp=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1); %Add graphics object to animate

        for qp=1:1:numel(hp) %For all graphics objects e.g. triangles/quads
            %         hp(qp).Marker=".";
            %         hp(qp).MarkerSize=markerSize2;
            hp(qp).FaceColor='interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(E_stress_mat_VM(:))]);
        axis(axisLim(V_DEF)); %Set axis limits statically
        view(140,30);
        camlight headlight;

        %     % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments

            % Compute 90% stress (area-weighted)
            sigmaVM = E_stress_mat_VM(:,:,qt); %von Misses stress for each elements
            Q = (sigmaVM(:) .* A) ./ A_tot; 
            p90_vec(qt) = prctile(Q,90);

            [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,qt));

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp(1) hp(1) hp(2) hp(2)]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        figure; hold on;
        plot(lambda_y_vec, p90_vec, 'r-o', 'LineWidth', 3, 'MarkerFaceColor', 'r');
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('90% Stress (Normalized)', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('p_{90} Stress vs Axial Stretch  |  Fiber Angle = %d°', ang), 'FontSize', 20, 'FontWeight', 'bold');
        grid on; axis square;
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');

        hf=cFigure;
        tiledlayout(1,2);

        hf = cFigure;
        tiledlayout(1,2);

        % === Left: Max Stress Animation ===
        ax1 = nexttile;
        title(ax1, ...
    sprintf('Max Stress $\\sigma_{vm}$ [MPa] | Fiber Angle = %d° | Elements = %d', ...
    ang, size(E,1)), ...
    'Interpreter','latex','FontSize',18);
        hp1 = gpatch(Fb, V_DEF(:,:,end), CV, 'none', 1, lineWidth1);
        axisGeom; colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(E_stress_mat_VM(:))]);
        view(140,30); camlight headlight;

        % === Right: 90% Stress Animation ===
        ax2 = nexttile;
      title(ax2, ...
    sprintf('90%% Stress $\\sigma_{vm}$ [MPa] | Fiber Angle = %d° | Elements = %d', ...
    ang, size(E,1)), ...
    'Interpreter','latex','FontSize',18);
        hp2 = gpatch(Fb, V_DEF(:,:,end), CV, 'none', 1, lineWidth1);
        axisGeom; colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(p90_vec)]);
        view(140,30); camlight headlight;




        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        CV=abs(N_disp_mat(:,1,end)); %Current displacement magnitude

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy

        gtitle([febioFebFileNamePart,': Press play to animate']);
        title(sprintf('$u_x$ [mm]  |  Fiber Angle = %d°', ang), 'Interpreter', 'latex');
        hp=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1); %Add graphics object to animate

        for qp=1:1:numel(hp) %For all graphics objects e.g. triangles/quads
            %         hp(qp).Marker=".";
            %         hp(qp).MarkerSize=markerSize2;
            hp(qp).FaceColor='interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0     2.0668]);
        % caxis([0 max(CV)]);
        % axis(axisLim(V_DEF)); %Set axis limits statically
        axis([-5 30 -5 90 -10 10])
        view(140,30);
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments

            CV=abs(N_disp_mat(:,1,qt)); %Current displacement magnitude

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp(1) hp(1) hp(2) hp(2)]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;
    end

end
%%
function [FT,VT,CT]=replicatemesh(F,V,p,d)
numcopies=prod(p);
numfaces=length(F);
numvertices=length(V);
C=ones(numfaces,1);
FT=repmat(F,numcopies,1);
VT=repmat(V,numcopies,1);
CT=repmat(C,numcopies,1);
c=1;
cf=1;
cv=1;
indexoffset=0;
for i=1:1:p(1)
    for j=1:1:p(2)
        FT(cf:cf+numfaces-1,:)=FT(cf:cf+numfaces-1,:)+indexoffset;
        VT(cv:cv+numvertices-1,:)=VT(cv:cv+numvertices-1,:)+[(i-1)*d(1),(j-1)*d(2)];
        CT(cf:cf+numfaces-1,:)=c;
        indexoffset=indexoffset+numvertices;
        c=c+1;
        cf=cf+numfaces;
        cv=cv+numvertices;
    end
end

[FT,VT]=mergeVertices(FT,VT);

end


%% plot
% clear; close all; clc; gohFolder='C:\Sahar-GIBBON-codes\paper3-final\final results-GOH\function-2fibers\data'; ogdenFolder='C:\Sahar-GIBBON-codes\paper3-final\final results-GOH\ogden-final\data'; geomIDs=[1 3 4 5]; geomNames={'Alternative slits (AS)-1st order','I re-entrant slits (IS)','Alternative slits (AS)-2nd order','H re-entrant (HS)'}; geomShort={'AS-1','IS','AS-2','HS'}; angles=0:15:90; fontSize=16; numPlotPointsGraphs=100; colors=[0 0 0;0.8500 0.3250 0.0980;0 0.4470 0.7410;0.4660 0.6740 0.1880;0.4940 0.1840 0.5560;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0.9290 0.6940 0.1250]; markers={'o','s','^','d','v','>','<','p'}; legendNames={'Ogden','GOH (\theta = 0°)','GOH (\theta = 15°)','GOH (\theta = 30°)','GOH (\theta = 45°)','GOH (\theta = 60°)','GOH (\theta = 75°)','GOH (\theta = 90°)'}; FendGOH=zeros(4,7); FendOgden=zeros(4,1); for g=1:4, figure; hold on; Hn=gobjects(8,1); D=dir(fullfile(ogdenFolder,sprintf('geometry_%d',geomIDs(g)),'*force_out.txt')); A=readmatrix(fullfile(D(1).folder,D(1).name),'FileType','text'); A=A(all(~isnan(A(:,1:4)),2),1:4); id=A(:,1); s=[1;find(diff(id)<=0)+1]; e=[s(2:end)-1;size(A,1)]; Fy=arrayfun(@(q)abs(sum(A(s(q):e(q),3))),1:numel(s)); lambdaRaw=linspace(1,1.25,numel(Fy)); lambdaPlot=linspace(1,1.25,numPlotPointsGraphs); FyPlot=interp1(lambdaRaw,Fy,lambdaPlot,'pchip','extrap'); FendOgden(g)=Fy(end); Hn(1)=plot(lambdaPlot,FyPlot,'Color',colors(1,:),'LineWidth',1.5,'Marker',markers{1},'MarkerIndices',1:10:length(lambdaPlot),'MarkerSize',6); for k=1:7, fileName=fullfile(gohFolder,sprintf('geometry_%d',geomIDs(g)),sprintf('geom%d_ang_%03d_force_out.txt',geomIDs(g),angles(k))); A=readmatrix(fileName,'FileType','text'); A=A(all(~isnan(A(:,1:4)),2),1:4); id=A(:,1); s=[1;find(diff(id)<=0)+1]; e=[s(2:end)-1;size(A,1)]; Fy=arrayfun(@(q)abs(sum(A(s(q):e(q),3))),1:numel(s)); lambdaRaw=linspace(1,1.25,numel(Fy)); lambdaPlot=linspace(1,1.25,numPlotPointsGraphs); FyPlot=interp1(lambdaRaw,Fy,lambdaPlot,'pchip','extrap'); FendGOH(g,k)=Fy(end); Hn(k+1)=plot(lambdaPlot,FyPlot,'Color',colors(k+1,:),'LineWidth',1.5,'Marker',markers{k+1},'MarkerIndices',1:10:length(lambdaPlot),'MarkerSize',6); end; title([geomNames{g} '|Force vs. Axial Stretch'],'FontSize',fontSize,'FontWeight','bold'); xlabel('\lambda_y (Axial Stretch)','FontSize',fontSize,'FontWeight','bold'); ylabel('|F_y|','FontSize',fontSize,'FontWeight','bold'); legend(Hn,legendNames,'Interpreter','tex','Location','northwest'); axis tight; axis square; box on; grid on; set(gca,'FontSize',fontSize,'FontWeight','bold','XMinorTick','on','YMinorTick','on'); drawnow; end
% theta_full=0:15:360; geomColors=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880;0.4940 0.1840 0.5560]; figure; hold on; for g=1:4, f=FendGOH(g,:); f360=[f fliplr(f(2:end-1)) f fliplr(f(2:end-1)) f(1)]; plot(theta_full,f360,'LineWidth',3,'Color',geomColors(g,:),'Marker','o','MarkerIndices',1:2:numel(theta_full),'MarkerSize',6,'DisplayName',[geomNames{g} ' GOH']); plot(theta_full,FendOgden(g)*ones(size(theta_full)),'--','LineWidth',2,'Color',geomColors(g,:),'DisplayName',[geomNames{g} ' Ogden']); end; xlabel('Fiber Angle (deg)','FontSize',16,'FontWeight','bold'); ylabel('|F_y| at \lambda_y = 1.25','FontSize',16,'FontWeight','bold'); title('Force vs Orientation (0–360°)','FontSize',22,'FontWeight','bold'); xlim([0 360]); xticks(0:30:360); legend('Location','eastoutside'); axis square; box on; grid on; set(gca,'YMinorTick','on','XMinorTick','on','FontSize',16,'FontWeight','bold')
% 
% theta_full=0:15:360; theta_circle=linspace(0,360,361); for g=1:4, f=FendGOH(g,:); f360=[f fliplr(f(2:end-1)) f fliplr(f(2:end-1)) f(1)]; figure; pax=polaraxes; pax.Layer='top'; hold on; polarplot(deg2rad(theta_full),f360,'LineWidth',2,'Color',colors(2,:)); polarplot(deg2rad(theta_circle),FendOgden(g)*ones(size(theta_circle)),'LineWidth',2,'Color',colors(1,:)); pax.ThetaTick=0:30:330; pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise'; pax.FontSize=16; pax.FontWeight='bold'; title([geomNames{g} ' | Polar Force at \lambda_y = 1.25'],'FontSize',22,'FontWeight','bold'); legend('GOH','Ogden','Location','best'); rmin=min([f360(:); FendOgden(g)]); rmax=max([f360(:); FendOgden(g)]); offset=0.01*(rmax-rmin); rlim([rmin-offset rmax+offset]); grid on; end
