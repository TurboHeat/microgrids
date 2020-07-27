%--------------------------------------------------------------%
% File: generate_all_figures.m (script)
% Author: Miguel Dias
% Date 25/08/16
% v1.0
% Description: Analysis of the solution found by the shortest_path solver:
% Plots the figures of tariff, demand and solution for all 36 cases.
% Exports the figures to eps format so they can be compiled into a single
% document using LaTex.
%--------------------------------------------------------------%
close all; clearvars; clc;
addpath(genpath('C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data'));
load graph_data_all_days.mat ;

% d_index; 1-3-> Large; 4-6->FSR; 7-9->Small restaurant; % 10-12->Residential
%d_index=mod(index,12)+[mod(index,12)==0]*12;
end_time = 24; %end time in h
dt = 15; %in s
n_lines = 3600 / dt;
T = end_time * n_lines; %number of time steps
t_plot = (1:n_lines * end_time)' ./ n_lines;
savepath = 'C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data\Case_Study_Images\';
save_fig = 1;

%%
joule2kWh = 1 / 3.6e6; %1kWh=3.6e6J
fuel_index = [1 * ones(12, 1); 2 * ones(12, 1); 3 * ones(12, 1)];

%% Initialize all vectors for the loop
electric_energy_kWh = zeros(1, numel(fuel_index));
heat_energy_kWh = zeros(1, numel(fuel_index));
base_electricity_charge = zeros(1, numel(fuel_index));
base_heat_charge = zeros(1, numel(fuel_index));
total_charge = zeros(1, numel(fuel_index));
power_MGT = zeros(T, numel(fuel_index));
heat_MGT = zeros(T, numel(fuel_index));
mdot_MGT = zeros(T, numel(fuel_index));
new_demand = zeros(T, numel(fuel_index));
path_cost = zeros(1, numel(fuel_index));
bought_elec = zeros(1, numel(fuel_index));
sold_energy = zeros(1, numel(fuel_index));
bought_fuel = zeros(1, numel(fuel_index));
bought_heat = zeros(1, numel(fuel_index));
MGT_cost = zeros(1, numel(fuel_index));
savings = zeros(1, numel(fuel_index));
FC = zeros(1, numel(fuel_index));
MGT_PDC = zeros(1, numel(fuel_index));
ut_PDC = zeros(1, numel(fuel_index));
MGT_IDC = zeros(1, numel(fuel_index));
ut_IDC = zeros(1, numel(fuel_index));

%%
FS = 20;
for i = 1:numel(fuel_index)
  d_index = mod(i, 12) + (mod(i, 12) == 0) * 12;
  %Compute base costs
  electric_energy_kWh(:, i) = sum(power_demand(:, d_index).*dt) * joule2kWh; %sum(w*s)= energy in kWh
  heat_energy_kWh(:, i) = sum(heat_demand(:, d_index)*dt) * joule2kWh;
  base_electricity_charge(:, i) = sum((power_demand(:, d_index) .* dt).*joule2kWh.*elec_tariff(:, d_index));
  base_heat_charge(:, i) = sum((heat_demand(:, d_index) .* dt).*joule2kWh.*heat_tariff(fuel_index(i)));
  total_charge(:, i) = base_electricity_charge(:, i) + base_heat_charge(:, i);
  
  figure(1);
  plot(t_plot, power_demand(:, d_index)/1e3, t_plot, heat_demand(:, d_index)/1e3, 'r'); % kWh
  legend('Electricity', 'Heat');
  xlabel('Time (h)');
  ylabel('Demanded power(kW)');
  title('Electricity/Heat energy demand')
  figure(2);
  plot(t_plot, elec_tariff(:, d_index));
  xlabel('Time (h)');
  ylabel('Rate ($/kWh)');
  title('Applicable electric rate');
  
  g = digraph(state_from, state_to, decided_costs(:, i));
  [path, path_length] = shortestpath(g, 'Start', 'End', 'Method', 'acyclic');
  [power_MGT(:, i), heat_MGT(:, i), mdot_MGT(:, i)] = extract_path(path, power_map, heat_map, fuel_map, SV_states);
  path_cost(:, i) = path_length;
  
  h1 = figure(3);
  plot(t_plot, power_demand(:, d_index)/1e3, t_plot, power_MGT(:, i)/1e3, t_plot, (power_demand(:, d_index) - power_MGT(:, i))/1e3, 'LineWidth', 2); %in kW
  set(gca, 'fontsize', 14)
  %legend('Power demand','CHP commitement', 'Utility commitement');
  xlim([0, 25]);
  ylim([0, Inf]);
  grid on;
  xlabel('Time [h]', 'FontSize', FS);
  ylabel('Demand and commitment [kW]', 'FontSize', FS);
  hold off;
  
  h2 = figure(4);
  plot(t_plot, heat_demand(:, d_index)/1e3, t_plot, heat_MGT(:, i)/1e3, t_plot, (heat_demand(:, d_index) - heat_MGT(:, i))/1e3, 'LineWidth', 2); %in Kw
  set(gca, 'fontsize', 14)
  %legend('Heat demand','CHP commitement', 'Utility commitement');
  xlim([0, 25]);
  ylim([0, Inf]);
  grid on;
  xlabel('Time [h]', 'FontSize', FS);
  ylabel('Demand and commitment [kW]', 'FontSize', FS);
  
  %{
    if save_fig
        figname=['FC', num2str(fuel_index(i)), 'B', num2str(tariff_map(d_index,1)), 'D', num2str(tariff_map(d_index,2))];
        print('-f1',[savepath, 'PH_Demand', figname],'-depsc');
        print('-f2',[savepath, 'Tariff', figname],'-depsc' );
        print('-f3',[savepath, 'Power_MGT', figname],'-depsc');
        print('-f4',[savepath, 'Heat_MGT', figname],'-depsc');
    end
  %}
  %Save just case study .fig files and eps
  if save_fig && fuel_index(i) == 1
    figname = ['FC', num2str(fuel_index(i)), 'B', num2str(tariff_map(d_index, 1)), 'D', num2str(tariff_map(d_index, 2))];
    savefig(h1, [savepath, 'Power_MGT', figname, '.fig']);
    savefig(h2, [savepath, 'Heat_MGT', figname, '.fig']);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-f3', [savepath, 'eps_noleg\', 'Power_MGT', figname], '-depsc')
    print('-f4', [savepath, 'eps_noleg\', 'Heat_MGT', figname], '-depsc')
  end
  new_demand(:, i) = subplus(power_demand(:, d_index)-power_MGT(:, i));
  bought_elec(i) = sum(subplus(power_demand(:, d_index)-power_MGT(:, i)).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
  sold_energy(i) = sum(subplus(-1.*(power_demand(:, d_index) - power_MGT(:, i))).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
  bought_fuel(i) = sum(mdot_MGT(:, i)*dt*price_kg_f(fuel_index(i))); %in $
  bought_heat(i) = sum(subplus(heat_demand(:, d_index)-heat_MGT(:, i)).*dt*joule2kWh.*heat_tariff(fuel_index(i))); %in $
  MGT_cost(i) = bought_elec(i) - sold_energy(i) + bought_fuel(i) + bought_heat(i); %in $
  [MGT_PDC(i), MGT_IDC(i), ut_PDC(i), ut_IDC(i), FC(i)] = GenerateDemandCharges(d_index, dt, power_demand(:, d_index), new_demand(:, i));
  % Absolute savings:
  savings(:, i) = total_charge(:, i) - path_cost(:, i);
end
disp(['All data saved to folder ', savepath]);

%% Economic metrics:
%t_path=[path.', circshift(path,1,2).'];

%% Save economic data
savepath = 'C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data\';
save([savepath, 'econ_data.mat'], 'elec_tariff', 'heat_tariff', ...
  'power_demand', 'heat_demand', 'power_MGT', 'heat_MGT', 'MGT_cost', ...
  'electric_energy_kWh', 'heat_energy_kWh', 'base_electricity_charge', ...
  'base_heat_charge', 'total_charge', 'bought_elec', 'sold_energy', ...
  'bought_heat', 'bought_fuel', 'MGT_cost', 'savings', 'new_demand', ...
  'FC', 'ut_PDC', 'ut_IDC', 'MGT_IDC', 'MGT_PDC', 'mdot_MGT');