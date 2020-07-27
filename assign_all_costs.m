%--------------------------------------------------------------%
% File: assign_all_costs.m (script)
% Author: Miguel Dias
% Date 14/08/16
% v1.0
% Description: Run all combinations of days and building to check if the
% tariffs are beeing correcty applied to all days of operations.
% Plots all combinations. The tariffs are valid for all seasons.
% Saves the tariffs and power demand combinations in a mat file, along with
% a map of the shape (BUILDING, DAY) (columns) and the corresponding linear
%index in the row.
% Also runs Iliya's mapping to save all needed variables in one file
%--------------------------------------------------------------%
clearvars;
close all;
clc;
addpath(genpath('C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data'));

%%
show_plot = 0; %boolean to decide if the plots off tariffs and demands are shown
smooth_demand = 1; % decide to perform or not smoothing on the demand data
smooth_time = 21; %number of steps for smoothing, needs to be odd. 15=3.75min, 21=5.25min
dt = 15;
n_lines = 3600 / dt; %number of time steps in 1h
end_time = 24; %in hours
T = end_time * n_lines; %number of time steps
BUILDING = [1, 2, 3, 4];
DAY      = [1 2 3];
total_nodes = 46; %smax*vmax+1
buildingstring = {'Large Hotel', 'Full Service Restaraunt', 'Small Hotel', 'Residential'}; %for plotting
daystring = {'Winter', 'Spring/Autumn', 'Summer'}; %for plotting
coststring = {'C1', 'C2', 'C3'}; %{C1, C2 C3}=[7.74 8.85 6.80] $/1000ft^3

%% Define save folder
cname = getenv('computername');
if strcmp(cname, 'AEJETLAB38') %at uni
  savepath = 'C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data\';
else %my PC
  savepath = 'C:\Users\migue\Dropbox\Technion Grand Energy Program\Miguel Dias\Data\';
end

%% get all tariffs combinations

%[Building; day]
tariff_map = uint8([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4; 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]).';
elec_tariff = zeros(T, numel(BUILDING)*numel(DAY));
power_demand = zeros(T, numel(BUILDING)*numel(DAY));
heat_demand = zeros(T, numel(BUILDING)*numel(DAY));
i = 1;
for b = 1:length(BUILDING)
  for d = 1:length(DAY)
    elec_tariff(:, i) = CreateElectricityTariffProfile(b, d, dt);
    [power_demand(:, i), heat_demand(:, i)] = CreateDemandProfileVector(b, d, dt);
    % SMOOTHING
    if smooth_demand
      power_demand(:, i) = smooth(power_demand(:, i), smooth_time); %select appropriate amount of hours
      heat_demand(:, i) = smooth(heat_demand(:, i), smooth_time);
    end
    i = i + 1;
  end
end
heat_tariff = 0.04; %define heat tariff

%% plotting according to building type
if show_plot
  t = [1:n_lines * end_time]' ./ n_lines;
  %Plot tariffs
  k = 1;
  for i = 1:4 %building
    figure(i);
    suptitle([buildingstring{i}; {''}; {''}]); %empty chars to get correct spacing
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
  for i = 1:4 %building
    figure(i+4);
    suptitle([buildingstring{i}; {''}; {''}]); %empty chars to get correct spacing
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
load graph24h.mat 'Transition_t';
state_from = Transition_t.state_t;
state_to = Transition_t.state_t1;
[SV_states, time_from, n_tsteps, from_state_map, to_state_map] = new_mappings_iliya(Transition_t);
clear Transition_t;
load 'sv_mappings'

%% Natural gas parameter
Qr = 49736500; %J/kg
h_env = 3.015e5; %J/kg
h_100 = 3.9748e5; %J/kg
eta_HRU = 0.89;
eta_b = .98;
mair_mf = 17.2 * 1.2;
Ph_mf = eta_HRU * (mair_mf + 1) * ((Qr * eta_b + mair_mf * h_env) / (mair_mf + 1) - h_100) / 3.6e6; %kWh/kg
price_ft3 = [7.74, 8.85, 6.80]; % $/1000ft^3
density_CH4 = 0.68; % kg/m^3
ft3_m3 = power(0.3048, 3);
price_m3 = price_ft3 ./ (1000 * ft3_m3); %price in $/m^3
price_kg_f = price_m3 / density_CH4; % for MGT costs
price_kWh = price_kg_f / Ph_mf; % price in $/kWh, for heat tariff
heat_tariff = price_kWh;

%% Reshape power maps
power_map = [reshape(power_map.', [], 1); 0];
heat_map = [reshape(heat_map.', [], 1); 0];
fuel_map = [reshape(m_dot_f_map.', [], 1); 0];
mdot_fuel_SU = sum(m_dot_f_map(:, 1:end-1)+m_dot_f_map(:, 2:end), 2); % Computes mdot of fuel consumed in startup for the 5 different max states
% Select right column with the final price
sol_select = [~SV_states(from_state_map, 1) & ~SV_states(to_state_map, 1), ... % Off-off
  ~SV_states(from_state_map, 1) & SV_states(to_state_map, 1), ... % Start up
  ~SV_states(to_state_map, 1), ... % Shut down
  true(numel(n_tsteps), 1)];% Remaining transitions
[~, sol_select] = max(sol_select, [], 2);
% assigns a small penalty to every input (s,v) change
transition_penalty = [zeros(total_nodes, 1); ...
  ~[SV_states(from_state_map(total_nodes+1:end-total_nodes), 1) ==, SV_states(to_state_map(total_nodes+1:end-total_nodes), 1) & ... %checks equality of S values
  SV_states(from_state_map(total_nodes+1:end-total_nodes), 2) ==, SV_states(to_state_map(total_nodes+1:end-total_nodes), 2)]; ...
  zeros(total_nodes, 1)];%checks equality of V values
%% Main loop to assign edges
k = 1;
decided_costs = zeros(length(to_state_map), numel(BUILDING)*numel(DAY)*numel(price_kg_f));
for cost = 1:numel(price_kg_f)
  fuel_price = price_kg_f(cost);
  for i = 1:numel(BUILDING) * numel(DAY)
    decided_costs(:, k) = assign_costs_fcn(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, total_nodes, sol_select, time_from, n_tsteps, from_state_map, to_state_map, power_demand(:, i), heat_demand(:, i), elec_tariff(:, i), heat_tariff(cost), fuel_price, transition_penalty);
    k = k + 1;
  end
end

%%
save([savepath, 'all_data.mat'], 'tariff_map', 'elec_tariff', 'power_demand', ...
  'heat_demand', 'SV_states', 'time_from', 'n_tsteps', 'heat_tariff', ...
  'from_state_map', 'to_state_map', 'heat_map', 'm_dot_f_map', 'power_map', ...
  'price_kg_f', 'decided_costs', 'state_from', 'state_to');

%g1=digraph(state_from, state_to, decided_costs(:,1));
save([savepath, 'graph_data_all_days.mat']);