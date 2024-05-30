run("startup.m")

%% Add necessary MRST modules
mrstModule add ad-props ad-core ad-blackoil compositional geothermal mrst-gui

%% Set up the grid
cartDim = [60, 60, 60];
L       = [600, 600, 600];
G = cartGrid(cartDim, L);
G = computeGeometry(G);

%% Load data
data = importdata('MRST_Perm_Poro_BestFit.csv').data;
poro = data(:,4);
perm_i = data(:,5);
perm_j = data(:,6);
perm_k = data(:,7);

rock.poro = poro;
rock.perm = [perm_i perm_j perm_k];


%% Plot porosity
figure; 
plotCellData(G, rock.poro); 
% view([0 180]);
colorbar;
set(gcf, 'position', [100,300,550,400]);
title('Porosity');
axis tight;
set(gca,'XDir','normal')  
set(gca,'YDir','normal') 
set(gca,'ZDir','normal')

%% Set fluid structure properties
rhoWS = 1000;
% Define fluid structure
fluid = initSimpleADIFluid('mu'    , 1.0e-3, ...
                           'rho'   , rhoWS , ...
                           'phases', 'W'   );
fluid = addThermalFluidProps(fluid           , ... % Original fluid
                             'Cp'     , 4.2e3, ... % Specific heat capacity
                             'lambdaF', 0.6  , ... % Thermal conductivity
                             'useEOS' , true );    % Use equation of state

%% Make rock structure
% add thermal props
rock = addThermalRockProps(rock           , ... % Original rock
                           'CpR'    , 790, ... % Specific heat capacity
                           'lambdaR', 3.05   , ... % Thermal conductivity
                           'rhoR'   , 2640, ... % Rock density
                           'tau'    , 1   );    % Tortuosity

%% Define the model
gravity reset on
gravity([0, 0, 9.81]);  
model = GeothermalModel(G, rock, fluid);

%% Set up initial state
state0   = initResSol(G, 1, 1);
initial_data = importdata('MRST_BC_IC.csv').data;  % load data
state0.pressure = initial_data(:,4);  % pressure
% state0.pressure = 1.5*initial_data(:,4);  % pressure
state0.T = initial_data(:,5);         % temperature

%% Set Boundary conditions and schedule
bc = [];

initial_pressure = reshape(state0.pressure,cartDim(1,1),cartDim(1,2),cartDim(1,3));
initial_temperature = reshape(state0.T,cartDim(1,1),cartDim(1,2),cartDim(1,3));

pressure_zmin = initial_pressure(:,:,1);
pressure_zmax = initial_pressure(:,:,end);
temperature_zmin = initial_temperature(:,:,1);
temperature_zmax = initial_temperature(:,:,end);

pressure_zmin = reshape(pressure_zmin,cartDim(1,1)*cartDim(1,2),1);
pressure_zmax = reshape(pressure_zmax,cartDim(1,1)*cartDim(1,2),1);
temperature_zmin = reshape(temperature_zmin,cartDim(1,1)*cartDim(1,2),1);
temperature_zmax = reshape(temperature_zmax,cartDim(1,1)*cartDim(1,2),1);

bc  = pside(bc, G, 'ZMin', pressure_zmin ,'sat', 1); 
bc  = pside(bc, G, 'ZMax', pressure_zmax ,'sat', 1); 

Tbc = [temperature_zmin;temperature_zmax];
bc  = addThermalBCProps(bc, 'T', Tbc);          % Temperature BCs

%% Wells
nx = G.cartDims(1);
ny = G.cartDims(2);

% injection rate
inj_rate1 = 0.0065;   % 2.5 bpm
inj_rate2 = 0.013;   % 5 bpm

W = [];
% Injection Well
Inj_point = 45+26*ny+19*nx*ny;
W = addWell(W, G, rock, Inj_point, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate1, ... % volumetric injection rate
            'comp_i', [1 0], ...  % inject water
            'Sign',1);    

Tinj = (273.15 + 20)*Kelvin;             % 20 degrees Celsius
W    = addThermalWellProps(W, G, rock, fluid, 'T', Tinj); % Add temperature field T to wells

%% Schedule
schedule.control = struct('W', W, 'bc', bc);
schedule.control.W.val = 0;
schedule.control(2) = struct('W', W, 'bc', bc);
schedule.control(3) = struct('W', W, 'bc', bc);
schedule.control(4) = struct('W', W, 'bc', bc);
schedule.control(4).W.val = inj_rate2;

timesteps1 = 1*minute*ones(1, 1);
timesteps2 = 10*minute*ones(6, 1);
timesteps3 = 1*minute*ones(1, 1);
timesteps4 = 10*minute*ones(10, 1);
schedule.step.val = [timesteps1;timesteps2;timesteps3;timesteps4];
schedule.step.control = [ones(numel(timesteps1), 1);ones(numel(timesteps2),1)*2;    ...
                         ones(numel(timesteps3),1)*3; ones(numel(timesteps4),1)*4];

%% Run simulation
tic
[~, states] = simulateScheduleAD(state0, model, schedule);
time_spent = toc;

hour = floor(time_spent/3600);              
minute = floor(mod(time_spent,3600)/60);  
second = round(time_spent-3600*hour-60*minute);

disp(['Time Elapsed: ',mat2str(hour),' hr ',mat2str(minute),' min ',mat2str(second),' sec ']);

%% Plot Results
figure;
plotToolbar(G,states);
% view([0 180]);
colorbar;
set(gcf, 'position', [800,250,550,450]);
axis tight;
set(gca,'XDir','normal')  
set(gca,'YDir','normal') 
set(gca,'ZDir','normal')

save('Forge_CirculationTest_BestFit.mat','states');

%% Pressure change plots
obs = importdata('CirculationTest_Measurements.csv').data;
obs_time = obs(:,1);
obs_pressure = obs(:,2);

simul_time = cumsum(schedule.step.val)/60-10;
simul_time(1) = 0;

obs_num = length(simul_time);

simul_pressure = [];
for ns = 1:obs_num
    temp1 = states{ns}.pressure;
    temp2 = states{ns}.T;
    simul_pressure = [simul_pressure;temp1(Inj_point)];     % pressure data
end

simul_pressure = simul_pressure * 0.000145038;   % pa -> psi
simul_pressure = simul_pressure - simul_pressure(1);  % pressure change

figure;
plot(obs_time,obs_pressure);
hold on
plot(simul_time,simul_pressure);
title('Pressure Change');
xlabel('Time (min)') 
ylabel('Injection Relative Pressure (psi)') 
xlim([0 150])
hold off
legend('Measured','Simulated','Location','southeast');

save('data_pressure.txt','simul_pressure','-ascii');
