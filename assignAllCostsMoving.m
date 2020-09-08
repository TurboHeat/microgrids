function [varargout] = assignAllCostsMoving(kwargs)
%% Handle inputs
arguments
  kwargs.showPlot (1,1) logical = false % will plots of tariffs and demands be shown?
  kwargs.smoothDemandTimesteps (1,1) double {mustBeInteger, mustBePositive} = 21; % number of timesteps for smoothing. 1=off
  kwargs.endTime (1,1) double {mustBePositive} = 24; % [h]
  kwargs.savePath (1,1) string = "../Data/";
  kwargs.transitionPenalty (1,1) double = 0.01;
  kwargs.demandStandardEnvelope (1,1) double {mustBeNonnegative} = 0; % α in: expectedDemand = μ ± α·σ
  kwargs.timeStepSize (1,1) double = 15; % duration of time step in [s]
end
% Unpack kwargs:
showPlot = kwargs.showPlot;
smoothPeriod = kwargs.smoothDemandTimesteps;
endTime = kwargs.endTime;
savePath = kwargs.savePath;
transitionPenalty = kwargs.transitionPenalty;
alpha = kwargs.demandStandardEnvelope;
dt = kwargs.timeStepSize; 

%% Constants
SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;

% Natural gas parameters
[PRICE_kg_f, HEAT_TARIFF] = NATURAL_GAS_PARAMS();

% Time-related definitions
N_LINES = SECONDS_PER_HOUR / dt; % number of time steps in 1h
T = endTime * N_LINES; % total number of time steps

% Buildings
BUILDING = BuildingType(1:4);
BLDG_DESCR = {'Large Hotel', 'Full Service Restaraunt', 'Small Hotel', 'Residential'}; %for plotting
NUM_BUILDINGS = numel(BUILDING);

% Demands
[CHP, NUM_WINDOWS] = LOAD_DEMAND_DATASETS();

%% Get all tariff & demand combinations
elecTariffs = cell(NUM_BUILDINGS, NUM_WINDOWS);
demands = cell(NUM_BUILDINGS, NUM_WINDOWS);

q = 1;
for b = 1:numel(CHP)
  for d = 1:NUM_WINDOWS
    elec_tariff(:, q) = createElectricityTariffProfile(b, d, dt);
    [power_demand(:, q), heat_demand(:, q)] = createDemandProfileVector(b, d, dt);
    % SMOOTHING
    if smoothDemand
      power_demand(:, q) = smooth(power_demand(:, q), smoothTime); %select appropriate amount of hours
      heat_demand(:, q) = smooth(heat_demand(:, q), smoothTime);
    end
    q = q + 1;
  end
end

%% plotting according to building type
if showPlot
  t = (1:N_LINES * endTime).' ./ N_LINES; % what is it?
  %Plot tariffs
  k = 1;
  for q = 1:4 %building
    figure(q);
    suptitle([BLDG_DESCR{q}; {''}; {''}]); %empty chars to get correct spacing
    for j = 1:3 %day
      subplot(1, 3, j);
      plot(t, elec_tariff(:, k));
      title(EXAMPLE_DAY_DESCR{j});
      xlabel('Time (h)');
      ylabel('Rate($/kWh)')
      k = k + 1;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
  end
  %Plot demands
  k = 1;
  for q = 1:4 %building
    figure(q+4);
    suptitle([BLDG_DESCR{q}; {''}; {''}]); %empty chars to get correct spacing
    for j = 1:3 %day
      subplot(1, 3, j);
      plot(t, power_demand(:, k), t, heat_demand(:, k));
      title(EXAMPLE_DAY_DESCR{j});
      legend('Power Demand', 'Heat Demand');
      xlabel('Time (h)');
      ylabel('Power(kW)');
      k = k + 1;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
  end
end

%% Run Iliya's mapping and save all the variables at once
load('graph_24h.mat', 'g', 'svToStateNumber');
state_from = g.Edges.EndNodes(:,1);
state_to = g.Edges.EndNodes(:, 2);

%%%
total_nodes = size(svToStateNumber, 1); %(smax-smin)*(vmax-vmin+1)+1
[SV_states, time_from, n_tsteps, from_state_map, to_state_map] = ...
  transformAdjacencyMatrix(state_from, state_to, svToStateNumber);
clear A g;
%% Load (and optionally rename) mappings
if isfile(fullfile(savePath, 'CHP.mat')) % Engine similar to Capstone C65 (incl. extrapolation)
  tmp = load(fullfile(savePath, 'CHP.mat'));
  power_map = tmp.PowerSurf;
  heat_map = tmp.PheatSurfExt;
  m_dot_f_map = tmp.FFSurfExt;
else
  tmp = load(fullfile(savePath, 'sv_mappings.mat')); % "Virtual" engine
  power_map = tmp.power_map;
  heat_map = tmp.heat_map;
  m_dot_f_map = tmp.m_dot_f_map;
end
clear tmp
% Verify that we have the right mappings:
assert( isequal( numel(power_map), numel(heat_map), numel(m_dot_f_map), total_nodes-1 ), ...
  'Mismatch between matrix sizes! Try using a different ''sv_mappings'' file.');
%% Reshape the power maps
power_map = [reshape(power_map.', [], 1); 0];
heat_map = [reshape(heat_map.', [], 1); 0];
fuel_map = [reshape(m_dot_f_map.', [], 1); 0];
mdot_fuel_SU = sum(m_dot_f_map(:, 1:end-1)+m_dot_f_map(:, 2:end), 2); % Computes mdot of fuel consumed in startup for all the different max states
% Select right column with the final price
sol_select = [~SV_states(from_state_map, 1) & ~SV_states(to_state_map, 1), ... % Off-off
  ~SV_states(from_state_map, 1) & SV_states(to_state_map, 1), ... % Start up
  ~SV_states(to_state_map, 1), ... % Shut down
  true(numel(n_tsteps), 1)];% Remaining transitions
[~, sol_select] = max(sol_select, [], 2);
% assigns a small penalty to every input (s,v) change
transition_penalty_indicator = [zeros(total_nodes, 1); ...
  ~(SV_states(from_state_map(total_nodes+1:end-total_nodes), 1) == SV_states(to_state_map(total_nodes+1:end-total_nodes), 1) & ... %checks equality of S values
    SV_states(from_state_map(total_nodes+1:end-total_nodes), 2) == SV_states(to_state_map(total_nodes+1:end-total_nodes), 2));
   zeros(total_nodes,1)]; %checks equality of V values
%% Main loop to assign edges
k = 1;
decided_costs = zeros(length(to_state_map), numel(BUILDING)*numel(DAY)*numel(PRICE_kg_f));
progressbar('Price levels [kg*f]','Building*day combinations');
nPrices = numel(PRICE_kg_f);
nConfigs = numel(BUILDING) * numel(DAY);
for cost = 1:nPrices % this takes ~4min
  fuel_price = PRICE_kg_f(cost);
  for q = 1:nConfigs
    decided_costs(:, k) = assignCosts(double(dt), power_map, heat_map, fuel_map, ...
      mdot_fuel_SU, total_nodes, sol_select, time_from, n_tsteps,...
      from_state_map, to_state_map,...
      power_demand.mean(:,q) + alpha * power_demand.std(:,q),...
      heat_demand.mean(:,q) + alpha * heat_demand.std(:,q), ...
      elec_tariff(:, q), HEAT_TARIFF(cost), fuel_price, transition_penalty_indicator,transitionPenalty);
    k = k + 1;
    progressbar([], q/nConfigs);
  end
  progressbar(cost/nPrices, []);
end
progressbar(1); % close progressbar

%% Save results to disk
% This requires ~3.4GB of space!
aux = decided_costs;
decided_costs = single(decided_costs);
save(fullfile(savePath, 'all_data.mat'), 'tariff_map', 'elec_tariff', 'power_demand', ...
  'heat_demand', 'SV_states', 'time_from', 'n_tsteps', 'heat_tariff', ...
  'from_state_map', 'to_state_map', 'heat_map', 'm_dot_f_map', 'power_map', ...
  'price_kg_f', 'decided_costs', 'state_from', 'state_to', 'sol_select');

%g1=digraph(state_from, state_to, decided_costs(:,1));
save(fullfile(savePath, 'graph_data_all_days.mat'), '-regexp', '[^aux]');
% Assign output:
if nargout
  varargout{1} = aux;
else
  assignin('base', 'decided_costs', aux);
end

end

function [price_kg_f, price_kWh] = NATURAL_GAS_PARAMS()
% This function contains some computations using constants. There is no need to perform
% them each time, so instead the end result is returned.
%{
Qr = 49736500;   % [J/kg]
h_env = 3.015e5; % [J/kg]
h_100 = 3.9748e5;% [J/kg]
eta_HRU = 0.89;
eta_b = 0.98;
mair_mf = 17.2 * 1.2;
Ph_mf = eta_HRU * (mair_mf + 1) * ((Qr * eta_b + mair_mf * h_env) / (mair_mf + 1) - h_100) / 3.6e6; %kWh/kg
price_ft3 = [7.74, 8.85, 6.80]; % $/1000ft^3
density_CH4 = 0.68; % kg/m^3
ft3_m3 = power(0.3048, 3);
price_m3 = price_ft3 ./ (1000 * ft3_m3); %price in $/m^3
price_kg_f = price_m3 / density_CH4; % for MGT costs
price_kWh = price_kg_f / Ph_mf; % price in $/kWh, for heat tariff
%}
price_kg_f = [0.401964000624002,0.459610000713491,0.353146667214886];
price_kWh = [0.0350691841359548,0.0400984857368475,0.0308101359333970];
end

function [chp, nWindows] = LOAD_DEMAND_DATASETS()
% The code below creates 2-week averaging windows for the 5 building types, where
% the first window is [02-Jan-2004 00:00:00, 16-Jan-2004 00:00:00] (because we don't have
% data from the end of 2003 and we don't want to use the end of 2004 as a substitute).
%{
chp = [CHP2004Provider("../Data/RefBldgLargeHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgFullServiceRestaurantNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgSmallHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/USA_NY_New.York-Central.Park.725033_TMY3_HIGH.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY)];     

startDay = datetime(2004, 1, 16); % skip the first two weeks
arrayfun(@(x)x.fastForward(startDay, 'last'), chp);

nextCnt = 0; % Should be equal to: 365-14+1 = 352
while (chp(1).hasNext)
  nextCnt = nextCnt + 1;
  [~] = chp(1).next();
end
%}
chp = struct2array(load('../Data/CHP2004.mat', 'chp'));
nWindows = struct2array(load('../Data/CHP2004.mat', 'nextCnt'));
end
