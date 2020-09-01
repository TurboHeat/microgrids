%--------------------------------------------------------------%
% File: assignAllCosts.m
% Author: Miel Sharf
% Date 05/08/20
%
% Description: Run all combinations of days and building to check if the
% tariffs are beeing correcty applied to all days of operations.
% Plots all combinations. The tariffs are valid for all seasons.
% Saves the tariffs and power demand combinations in a mat file, along with
% a map of the shape (BUILDING, DAY) (columns) and the corresponding linear
% index in the row.
% Also runs Iliya's mapping to save all needed variables in one file.
%--------------------------------------------------------------%
function [varargout] = assignAllCosts(kwargs)
%% Handle inputs
arguments
  kwargs.showPlot (1,1) logical = false % will plots of tariffs and demands be shown?
  kwargs.smoothDemand (1,1) logical = true; % will smoothing be performed on the demand data?
  kwargs.smoothTime (1,1) double {mustBeInteger, mustBePositive} = 21; % number of steps for smoothing. 15=3.75min, 21=5.25min
  kwargs.endTime (1,1) double {mustBePositive} = 24; % [h]
  kwargs.savePath (1,1) string = "../Data/"; 
end
% Unpack kwargs:
showPlot = kwargs.showPlot;
smoothDemand = kwargs.smoothDemand;
smoothTime = kwargs.smoothTime;
endTime = kwargs.endTime;
savePath = kwargs.savePath;

%% Constants
SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;

% Natural gas parameters
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
heat_tariff = price_kWh;

%
dt = 15; % duration of time step in [s]
n_lines = SECONDS_PER_HOUR / dt; % number of time steps in 1h
T = endTime * n_lines; % total number of time steps
BUILDING = [1, 2, 3, 4];
DAY      = [1, 2, 3];
buildingstring = {'Large Hotel', 'Full Service Restaraunt', 'Small Hotel', 'Residential'}; %for plotting
daystring = {'Winter', 'Spring/Autumn', 'Summer'}; %for plotting

%% Get all tariff combinations
%[Building; day]
tariff_map = uint8([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4; 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]).';
elec_tariff = zeros(T, numel(BUILDING)*numel(DAY));
power_demand = zeros(T, numel(BUILDING)*numel(DAY));
heat_demand = zeros(T, numel(BUILDING)*numel(DAY));
q = 1;
for b = 1:length(BUILDING)
  for d = 1:length(DAY)
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
  t = (1:n_lines * endTime).' ./ n_lines; % what is it?
  %Plot tariffs
  k = 1;
  for q = 1:4 %building
    figure(q);
    suptitle([buildingstring{q}; {''}; {''}]); %empty chars to get correct spacing
    for j = 1:3 %day
      subplot(1, 3, j);
      plot(t, elec_tariff(:, k));
      title(daystring{j});
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
    suptitle([buildingstring{q}; {''}; {''}]); %empty chars to get correct spacing
    for j = 1:3 %day
      subplot(1, 3, j);
      plot(t, power_demand(:, k), t, heat_demand(:, k));
      title(daystring{j});
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
transition_penalty = [zeros(total_nodes, 1); ...
  ~(SV_states(from_state_map(total_nodes+1:end-total_nodes), 1) == SV_states(to_state_map(total_nodes+1:end-total_nodes), 1) & ... %checks equality of S values
    SV_states(from_state_map(total_nodes+1:end-total_nodes), 2) == SV_states(to_state_map(total_nodes+1:end-total_nodes), 2));
   zeros(total_nodes,1)]; %checks equality of V values
%% Main loop to assign edges
k = 1;
decided_costs = zeros(length(to_state_map), numel(BUILDING)*numel(DAY)*numel(price_kg_f));
progressbar('Price levels [kg*f]','Building*day combinations');
nPrices = numel(price_kg_f);
nConfigs = numel(BUILDING) * numel(DAY);
for cost = 1:nPrices % this takes ~4min
  fuel_price = price_kg_f(cost);
  for q = 1:nConfigs
    decided_costs(:, k) = assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, total_nodes, sol_select, time_from, n_tsteps, from_state_map, to_state_map, power_demand(:, q), heat_demand(:, q), elec_tariff(:, q), heat_tariff(cost), fuel_price, transition_penalty);
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
