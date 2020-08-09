%--------------------------------------------------------------%
% File: generate_all_figures_multi_unit_case.m (script)
% Author: Miel Sharf
% Date 08/08/2020
% v2.0
% Description: Solves the economic dispatch problem using a mixture of
% shortest-path solvers and a gradient descent approach for the
% multipliers.
% Plots the figures of tariff, demand and solution for all 36 cases.
% Exports the figures to eps format so they can be compiled into a single
% document using LaTex.
%--------------------------------------------------------------%
%Startup:
% The variable decided_costs takes too much room to be saved normally. If
% it already exists in the workspace, don't delete it. If it does not
% exist, load it and convert it to double.
if (exist('decided_costs', 'var'))
  clearvars - except decided_costs;
  close all;
  clc;
  addpath(genpath('C:\GasTurbinesProject\OneDrive_2020-04-27\Energy Project\Data'));
  load('graph_data_all_days.mat', '-regexp', '[^decided_costs]');
else
  clearvars;
  close all;
  clc;
  addpath(genpath('C:\GasTurbinesProject\OneDrive_2020-04-27\Energy Project\Data'));
  load graph_data_all_days.mat;
  decided_costs = double(decided_costs);
end
% Make sure that key parameters are defined
if ~exist('end_time', 'var')
  end_time = 24; %end time in h.
end
if ~exist('dt', 'var')
  dt = 15; %step length in seconds
end
if ~exist('n_lines', 'var')
  n_lines = 3600 / dt;
end
if ~exist('T', 'var')
  T = end_time * n_lines; %number of time steps
end
if ~exist('joule2kWh', 'var')
  joule2kWh = 1 / 3.6e6; %1kWh=3.6e6J
end
if ~exist('fuel_index', 'var')
  DayBuildingCombinations = length(DAY) * length(BUILDING);
  fuel_index = [1 * ones(DayBuildingCombinations, 1); 2 * ones(DayBuildingCombinations, 1); 3 * ones(DayBuildingCombinations, 1)]; %Type of fuel in various iterations.
end

%% Define key parameters for the script
% Physical and Algorithmical Parameters
n_MGTs = 4; %number of Multi-Gas Turbines.
StepSizeParameter = 1e-1; %step size parameter
IterationLimit = 1500;
BufferSize = 1; %Find the best solution within BufferSize last iterations.

% Logic Parameters
t_plot = (1:n_lines * end_time)' ./ n_lines; %time vector for plotting graphs.
%If and where to save the graphs
save_fig = 1;
savepath = 'C:\GasTurbinesProject\OneDrive_2020-04-27\Energy Project\Data\Case_Study_Images\';

FontSize = 20; %Font Size

%% Initialize all vectors for the loop

power_MGTs = cell(n_MGTs, 1);
heat_MGTs = cell(n_MGTs, 1);
mdot_MGTs = cell(n_MGTs, 1);
graphs_MGTs = cell(n_MGTs, 1);

for mgt = 1:n_MGTs
  power_MGTs{mgt} = zeros(T, length(fuel_index));
  heat_MGTs{mgt} = zeros(T, length(fuel_index));
  mdot_MGTs{mgt} = zeros(T, length(fuel_index));
end

path_cost = zeros(n_MGTs, length(fuel_index));

power_MGT_Total = zeros(T, 1);
heat_MGT_Total = zeros(T, 1);
base_electricity_charge = zeros(1, length(fuel_index));
base_heat_charge = zeros(1, length(fuel_index));
total_charge = zeros(1, length(fuel_index));
new_demand = zeros(T, length(fuel_index));
bought_elec = zeros(1, length(fuel_index));
sold_energy = zeros(1, length(fuel_index));
bought_fuel = zeros(1, length(fuel_index));
bought_heat = zeros(1, length(fuel_index));
MGT_cost = zeros(1, length(fuel_index));
savings = zeros(1, length(fuel_index));
FC = zeros(1, length(fuel_index));
MGT_PDC = zeros(1, length(fuel_index));
ut_PDC = zeros(1, length(fuel_index));
MGT_IDC = zeros(1, length(fuel_index));
ut_IDC = zeros(1, length(fuel_index));

%% Trick for Graphs of Tubrines
% It is much quicker to change the edge weights by calling g.Edges.Weight
% directly rather than rebuilding the graph g by digraph(). However, the
% way that MATLAB saves the array of edge weights does not correspond to
% state_from and state_to. We find out the corresponding permutation by
% asking MATLAB to build the same digraph with weights 1,2,3,..., which
% serve as indices.

g = digraph(state_from, state_to, 1:size(decided_costs, 1));
EdgePermutationMap = g.Edges.Weight;
%Build a new graph with the same transitions and costs Costs by:
%  g_new = g;
%  g_new.Edges.Weight = Costs(EdgePermutationMap);

for mgt = 1:n_MGTs
  graphs_MGTs{mgt} = g;
  %Update the weights later.
end

%% Run the solution method - iterate over all day-building-fuel cost combinations.
for i = 1:length(fuel_index)
  d_index = mod(i, 12) + (mod(i, 12) == 0) * 12;
  %d_index is the day-building combination:
  %1 - Large Hotel, Winter
  %2 - Large Hotel, Transition Season
  %3 - Large Hotel, Summer
  %4 - Restaurant, Winter
  %5 - Restaurant, Transition Season
  %6 - Restaurant, Summer
  %7 - Small Hotel, Winter
  %8 - Small Hotel, Transition Season
  %9 - Small Hotel, Summer
  %10 - Residential, Winter
  %11 - Residential, Transition Season
  %12 - Residential, Summer
  
  StepSizeParameter = 1e-1; %step size parameter
  EndBuffer.Turbines = cell(n_MGTs, BufferSize);
  EndBuffer.TotalPower = zeros(T, BufferSize);
  EndBuffer.TotalHeat = zeros(T, BufferSize);
  EndBuffer.Price = zeros(1, BufferSize);
  
  
  lambda = zeros(T, 2); %Two multipliers for each time - one for power constraint and one for heat constraint.
  %     lambda(:,1) = elec_tariff(:,d_index); %power multipliers are bounded by power tariff from utility
  
  
  %Compute base costs
  base_electricity_charge(i) = sum((power_demand(:, d_index) * dt).*joule2kWh.*elec_tariff(:, d_index));
  base_heat_charge(i) = sum((heat_demand(:, d_index) * dt).*joule2kWh.*heat_tariff(fuel_index(i)));
  total_charge(i) = base_electricity_charge(i) + base_heat_charge(i);
  
  %Plot power demend, heat demend and electricity tariff as a function of
  %time.
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
  lambda_hist = zeros(IterationLimit, 1);
  
  for iter = 1:IterationLimit
    BufferIndex = mod(iter, BufferSize);
    if (BufferIndex == 0)
      BufferIndex = BufferSize;
    end
    decided_costs_extra = assign_costs_multi_unit_fcn_extra(total_nodes, sol_select, time_from, n_tsteps, from_state_map, to_state_map, power_map, heat_map, power_demand(:, d_index), heat_demand(:, d_index), lambda*joule2kWh*dt);
    overall_cost = decided_costs(:, i) + decided_costs_extra;
    
    %Every MGT solves the shortest path problem to minimize
    for mgt = 1:n_MGTs
      %ADD EFFECT OF LAMBDA!!!!!
      graphs_MGTs{mgt}.Edges.Weight = overall_cost(EdgePermutationMap) + 0.1 * is_transition(EdgePermutationMap);
      StartNode = 1;
      EndNode = max(state_to);
      [path, path_length] = shortestpath(graphs_MGTs{mgt}, StartNode, EndNode, 'Method', 'acyclic');
      [power_MGTs{mgt}(:, i), heat_MGTs{mgt}(:, i), mdot_MGTs{mgt}(:, i)] = extract_path(path, power_map, heat_map, fuel_map, SV_states);
      path_cost(mgt, i) = path_length;
    end
    for mgt = 1:n_MGTs
      EndBuffer.Turbines{mgt, BufferIndex}.Power = power_MGTs{mgt}(:, i);
      EndBuffer.Turbines{mgt, BufferIndex}.Heat = heat_MGTs{mgt}(:, i);
      EndBuffer.Turbines{mgt, BufferIndex}.mDot = mdot_MGTs{mgt}(:, i);
      EndBuffer.Turbines{mgt, BufferIndex}.path_cost = path_cost(mgt, i);
    end
    
    power_MGT_Total = 0;
    heat_MGT_Total = 0;
    for mgt = 1:n_MGTs
      power_MGT_Total = power_MGT_Total + power_MGTs{mgt}(:, i);
      heat_MGT_Total = heat_MGT_Total + heat_MGTs{mgt}(:, i);
    end
    EndBuffer.TotalPower(:, BufferIndex) = power_MGT_Total;
    EndBuffer.TotalHeat(:, BufferIndex) = heat_MGT_Total;
    
    %We need to minimize the cost of buying power and heat from the utility.
    %The cost function is of buying p_t units of power at time t from the utility
    %is powercost_t*p_t-lambda_{p,t}*p_t, where p_t >=0. If
    %powercost_t>lambda_{p,t}, the minimum is achieved at p_t = 0. Otherwise, the
    %minimum is minus infinity. Thus, the maximization problem on
    %lambda has constraints lambda_{p,t} <= powercost_t, and for these
    %values, we buy now power from the utility. A similar situation
    %happens for heat purchases...
    
    power_Utility = 0;
    heat_Utility = 0;
    
    
    %Compute gradients
    
    diff_lambda_p = power_demand(:, d_index) - power_MGT_Total - power_Utility;
    diff_lambda_h = heat_demand(:, d_index) - heat_MGT_Total - heat_Utility;
    
    % step-size for the subgradient algorithm
    
    delta_p = StepSizeParameter / sqrt(iter); %./(2*sqrt(iter).*abs(diff_lambda_p));
    delta_h = StepSizeParameter / sqrt(iter); %./(2*sqrt(iter).*abs(diff_lambda_h));
    
    % update
    lambda(:, 1) = lambda(:, 1) + delta_p * diff_lambda_p;
    lambda(:, 2) = lambda(:, 2) + delta_h * diff_lambda_h;
    
    %projected gradient descent - lambda <= cost of power and heat from
    %utility.
    if (sum(sum(lambda >= 0)) == 0 || sum(sum(lambda <= [elec_tariff(:, d_index), heat_tariff(floor((i - 1)/12)+1) * ones(size(lambda(:, 2)))])) == 0)
      StepSizeParameter = StepSizeParameter / 10;
    end
    
    
    lambda(:, 1) = min(lambda(:, 1), elec_tariff(:, d_index)); %power multipliers are bounded by power tariff from utility
    lambda(:, 2) = min(lambda(:, 2), heat_tariff(floor((i - 1)/12)+1)); %heat multipliers are bounded by heat tariff from utility
    
    
    %         lambda(:,1) = max(lambda(:,1),0); %power multipliers are bounded by power tariff from utility
    %         lambda(:,2) = max(lambda(:,2),0); %heat multipliers are bounded by heat tariff from utility
    if (sum(lambda(:, 1).*power_demand(:, d_index)+lambda(:, 2).*heat_demand(:, d_index)) < 0)
      lambda(:, 1) = 0;
      lambda(:, 2) = 0;
    end
    
    lambda_hist(iter) = lambda(1, 1);
    
    %Compute Price for buffer
    ElecBought = sum(subplus(power_demand(:, d_index)-power_MGT_Total).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
    EnergySold = sum(subplus(-1.*(power_demand(:, d_index) - power_MGT_Total)).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
    %         FuelBought = sum((mdot_MGTs{1}(:,i)+mdot_MGTs{2}(:,i)+mdot_MGTs{3}(:,i)+mdot_MGTs{4}(:,i))*dt*price_kg_f(fuel_index(i))); %in $
    FuelBought = sum((mdot_MGTs{1}(:, i) + mdot_MGTs{2}(:, i))*dt*price_kg_f(fuel_index(i))); %in $
    HeatBought = sum(subplus(heat_demand(:, d_index)-heat_MGT_Total).*dt*joule2kWh.*heat_tariff(fuel_index(i))); %in $
    EndBuffer.Price(BufferIndex) = ElecBought - EnergySold + FuelBought + HeatBought; %in $
  end
  
  %% Find best Solution from the Buffer - Deal with unstability due to insufficient iterations
  [BestPrice, BestIndex] = min(EndBuffer.Price);
  power_MGT_Total = EndBuffer.TotalPower(:, BestIndex);
  heat_MGT_Total = EndBuffer.TotalHeat(:, BestIndex);
  for mgt = 1:n_MGTs
    power_MGTs{mgt}(:, i) = EndBuffer.Turbines{mgt, BestIndex}.Power;
    heat_MGTs{mgt}(:, i) = EndBuffer.Turbines{mgt, BestIndex}.Heat;
    mdot_MGTs{mgt}(:, i) = EndBuffer.Turbines{mgt, BestIndex}.mDot;
  end
  toc;
  
  %% Plot
  if (0)
    h1 = figure(3);
    subplot(2, 1, 1);
    plot(t_plot, power_demand(:, d_index)/1e3, t_plot, power_MGT_Total/1e3, ...
      t_plot, (power_demand(:, d_index) - power_MGT_Total)/1e3, 'LineWidth', 2);
    legend('Power Demand', 'Turbines', 'Utility');
    ylabel('Power [kw]');
    xlabel('Time [sec]');
    subplot(2, 1, 2);
    plot(t_plot, power_MGTs{1}(:, i)/1e3, ...
      t_plot, power_MGTs{2}(:, i)/1e3, t_plot, power_MGTs{3}(:, i)/1e3, ...
      t_plot, power_MGTs{4}(:, i)/1e3, 'LineWidth', 2);
    
    h2 = figure(4);
    subplot(2, 1, 1);
    plot(t_plot, heat_demand(:, d_index)/1e3, t_plot, heat_MGT_Total/1e3, ...
      t_plot, (heat_demand(:, d_index) - heat_MGT_Total)/1e3, 'LineWidth', 2);
    legend('Heat Demand', 'Turbines', 'Utility');
    ylabel('Heat [kw]');
    xlabel('Time [sec]');
    subplot(2, 1, 2);
    plot(t_plot, heat_MGTs{1}(:, i)/1e3, ...
      t_plot, heat_MGTs{2}(:, i)/1e3, t_plot, heat_MGTs{3}(:, i)/1e3, ...
      t_plot, heat_MGTs{4}(:, i)/1e3, 'LineWidth', 2);
    
    h3 = figure(5);
    plot(1:length(lambda_hist), lambda_hist, 'LineWidth', 2);
    legend('\lambda_{t=1}^p vs. iteration');
    ylabel('\lambda');
    xlabel('Iteration');
    
    %%
    %Save just case study .fig files and eps
    if save_fig && fuel_index(i) == 1
      figname = ['FC', num2str(fuel_index(i)), 'B', num2str(tariff_map(d_index, 1)), 'D', num2str(tariff_map(d_index, 2))];
      savefig(h1, [savepath, 'Power_MGT', figname, '.fig']);
      savefig(h2, [savepath, 'Heat_MGT', figname, '.fig']);
      savefig(h3, [savepath, 'Multipliers', figname, '.fig']);
      set(gcf, 'PaperPositionMode', 'auto');
      print('-f3', [savepath, 'eps_noleg\', 'Power_MGT', figname], '-depsc')
      print('-f4', [savepath, 'eps_noleg\', 'Heat_MGT', figname], '-depsc')
      print('-f5', [savepath, 'eps_noleg\', 'Multipliers', figname], '-depsc')
    end
  end
  new_demand(:, i) = subplus(power_demand(:, d_index)-power_MGT_Total);
  bought_elec(i) = sum(subplus(power_demand(:, d_index)-power_MGT_Total).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
  sold_energy(i) = sum(subplus(-1.*(power_demand(:, d_index) - power_MGT_Total)).*dt*joule2kWh.*elec_tariff(:, d_index)); %in $
  %Make Generic for n_MGTs.
  bought_fuel(i) = sum((mdot_MGTs{1}(:, i) + mdot_MGTs{2}(:, i) + mdot_MGTs{3}(:, i) + mdot_MGTs{4}(:, i))*dt*price_kg_f(fuel_index(i))); %in $
  bought_heat(i) = sum(subplus(heat_demand(:, d_index)-heat_MGT_Total).*dt*joule2kWh.*heat_tariff(fuel_index(i))); %in $
  MGT_cost(i) = bought_elec(i) - sold_energy(i) + bought_fuel(i) + bought_heat(i); %in $
  [MGT_PDC(i), MGT_IDC(i), ut_PDC(i), ut_IDC(i), FC(i)] = GenerateDemandCharges(d_index, dt, power_demand(:, d_index), new_demand(:, i));
  % Absolute savings:
  %Miel: Consider the cost of bought power and heat from utility in
  %savings computation
  %     savings(:,i)=total_charge(:,i)-sum(path_cost(:,i));
  savings(i) = total_charge(i) - MGT_cost(i);
  
  disp(['Run Number #', int2str(i)]);
  disp(['Power Generation:  ', num2str(power_MGTs{1}(1, i))]);
  disp(['Heat Generation:  ', num2str(heat_MGTs{1}(1, i))]);
  disp(['--------------------------------']);
end
% toc
disp(['All data saved to folder ', savepath]);

%% Save economic data
savepath = 'C:\GasTurbinesProject\OneDrive_2020-04-27\Energy Project\Data\';
save([savepath, 'econ_data.mat'], 'elec_tariff', 'heat_tariff', ...
  'power_demand', 'heat_demand', 'power_MGTs', 'heat_MGTs', 'MGT_cost', ...
  'base_electricity_charge', ...
  'base_heat_charge', 'total_charge', 'bought_elec', 'sold_energy', ...
  'bought_heat', 'bought_fuel', 'MGT_cost', 'savings', 'new_demand', ...
  'FC', 'ut_PDC', 'ut_IDC', 'MGT_IDC', 'MGT_PDC', 'mdot_MGTs');