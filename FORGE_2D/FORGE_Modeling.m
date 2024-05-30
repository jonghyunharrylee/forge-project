run("startup.m")

%% Add necessary MRST modules
mrstModule add ad-props ad-core ad-blackoil compositional geothermal mrst-gui

%% Set up the grid
cartDim = [60, 1, 60];
L       = [600, 1, 600];
G = cartGrid(cartDim, L);
G = computeGeometry(G);

%% Load data
data = importdata('Upscaled_Fractures_Slice.csv').data;
poro = data(:,1);
perm_i = data(:,2);
perm_j = data(:,3);
perm_k = data(:,4);

rock.poro = poro;
rock.perm = [perm_i perm_j perm_k];

%% Plot porosity
figure; 
plotCellData(G, rock.poro); 
view([0 180]);
colorbar;
set(gcf, 'position', [100,300,550,400]);
title('Porosity')
axis tight;

%% Set fluid properties
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
initial_data = importdata('Initial.csv').data;  % load data
state0.pressure = initial_data(:,5);  % pressure
state0.T = initial_data(:,6);         % temperature

%% Set Boundary conditions
bc = [];

initial_pressure = reshape(state0.pressure,cartDim(1,1),cartDim(1,3));
initial_temperature = reshape(state0.T,cartDim(1,1),cartDim(1,3));

pressure_zmin = initial_pressure(1,:);
pressure_zmax = initial_pressure(end,:);
temperature_zmin = initial_temperature(1,:);
temperature_zmax = initial_temperature(end,:);

pressure_zmin = reshape(pressure_zmin,cartDim(1,1),1);
pressure_zmax = reshape(pressure_zmax,cartDim(1,1),1);
temperature_zmin = reshape(temperature_zmin,cartDim(1,1),1);
temperature_zmax = reshape(temperature_zmax,cartDim(1,1),1);

bc  = pside(bc, G, 'ZMin', pressure_zmin ,'sat', 1); 
bc  = pside(bc, G, 'ZMax', pressure_zmax ,'sat', 1); 

Tbc = [temperature_zmin;temperature_zmax];
bc  = addThermalBCProps(bc, 'T', Tbc);          % Temperature BCs

%% Wells
nz = G.cartDims(3);

% injection rate
inj_rate = 1e-4;

W = [];
% Well 1
W = addWell(W, G, rock, 20*nz+44, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [1 0]);    % inject water

% Well 2
W = addWell(W, G, rock, 24*nz+33, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [1 0]);    % inject water

% Well 3
W = addWell(W, G, rock, 30*nz+21, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [1 0]);    % inject water

Tinj = (273.15 + 20)*Kelvin;             % 20 degrees Celsius
W    = addThermalWellProps(W, G, rock, fluid, 'T', Tinj); % Add temperature field T to wells

%% Schedule
schedule.control = struct('W', W, 'bc', bc);

timesteps = rampupTimesteps(1500*day,30*day,2);
schedule.step.val = timesteps;
schedule.step.control = ones(numel(schedule.step.val), 1); 

%% Run simulation
[~, states] = simulateScheduleAD(state0, model, schedule);

%% Plot Results
figure;
plotToolbar(G,states);
view([0 180]);
colorbar;
set(gcf, 'position', [800,250,550,450]);
axis tight;
