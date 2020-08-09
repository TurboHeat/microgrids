%--------------------------------------------------------------%
% File: econ_analysis_tables.m (script)
% Author: Miguel Dias
% Date 30/08/16
% v1.0
% Description: Analysis of the solution found by the shorstest_path solver:
% 1) Savings;
% 2) Generates econ data tables, grouped by building and fuel costs and
% exports them to an excel spreadsheet
% The buildings are:
% 1. Large Hotel;
% 2. Full Service Restaurant;
% 3. Small Hotel;
% 4. Residential neighbourhood
%--------------------------------------------------------------%
close all; clearvars; clc;
addpath(genpath('C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data'));
load econ_data.mat ;
load graph_data_all_days.mat price_ft3;
fuel_index = [1 * ones(12, 1); 2 * ones(12, 1); 3 * ones(12, 1)];

canSave = 0; %Uses a try/catch mechanism - if we save too fast, we get an error...

%issue with fixed cost for res building:
FC([10:12, 22:24, 34:36]) = 1.65;

%% Total costs($/day of operation):
Total_utility_costs = total_charge + (ut_PDC + ut_IDC) / 30 + FC;
Total_MGT_costs = MGT_cost + (MGT_PDC + MGT_IDC) / 30 + FC;
abs_savings = Total_utility_costs - Total_MGT_costs;
%Absolute savings
RowNames = {'D1', 'D2', 'D3'};
LH = abs_savings([1:3; 13:15; 25:27]).'; %Rows: days, Columns: Fuel cost
FSR = abs_savings([4:6; 16:18; 28:30]).';
SH = abs_savings([7:9; 19:21; 31:33]).';
R = abs_savings([10:12; 22:24; 34:36]).';

abs_savings_t = table(LH, FSR, SH, R, 'RowNames', RowNames);
filename = 'econ_data.xlsx';

xlswrite(filename, {' Absolute savings all days, all fuel costs'}, 'Savings', 'A2');
writetable(abs_savings_t, filename, 'Sheet', 'Savings', 'Range', 'C3', 'WriteRowNames', true);

%% Save economic data into excel spreadsheet
RowNames = {'D1', 'D2', 'D3'};
LH = savings([1:3; 13:15; 25:27]).'; %Rows: days, Columns: Fuel cost
FSR = savings([4:6; 16:18; 28:30]).';
SH = savings([7:9; 19:21; 31:33]).';
R = savings([10:12; 22:24; 34:36]).';

savings_t = table(LH, FSR, SH, R, 'RowNames', RowNames);
filename = 'econ_data.xlsx';

xlswrite(filename, {' Fixed Pricing of energy savings all days, all fuel costs'}, 'Savings', 'A12');
writetable(savings_t, filename, 'Sheet', 'Savings', 'Range', 'C13', 'WriteRowNames', true);

%%
savings_PDC = ut_PDC - MGT_PDC;
savings_IDC = ut_IDC - MGT_IDC;

LH = [savings_PDC([1, 13, 25]); ...
  savings_IDC([1, 13, 25]); ...
  savings_PDC([2, 14, 26]); ...
  savings_IDC([2, 14, 26]); ...
  savings_PDC([3, 15, 27]); ...
  savings_IDC([3, 15, 27])];%Hospital columns: FC, rows: day {PDC IDC}

FSR = [savings_PDC([4, 16, 28]); ...
  savings_IDC([4, 16, 28]); ...
  savings_PDC([5, 17, 29]); ...
  savings_IDC([5, 17, 29]); ...
  savings_PDC([6, 18, 30]); ...
  savings_IDC([6, 18, 30])];%FSR columns: FC, rows: day {PDC IDC}

SH = [savings_PDC([7, 19, 31]); ...
  savings_IDC([7, 19, 31]); ...
  savings_PDC([8, 20, 32]); ...
  savings_IDC([8, 20, 32]); ...
  savings_PDC([9, 21, 33]); ...
  savings_IDC([9, 21, 33])];%SH columns: FC, rows: day {PDC IDC}

R = [savings_PDC([10, 22, 34]); ...
  savings_IDC([10, 22, 34]); ...
  savings_PDC([11, 23, 35]); ...
  savings_IDC([11, 23, 35]); ...
  savings_PDC([12, 24, 36]); ...
  savings_IDC([12, 24, 36])];%R columns: FC, rows: day {PDC IDC} NO PDC
%or IDC!!!
savings_d = table(LH, FSR, SH, R);
xlswrite(filename, {'Demand charges all days, all fuel costs'}, 'Savings', 'A18');
writetable(savings_d, filename, 'Sheet', 'Savings', 'Range', 'C19', 'WriteRowNames', true);

%
% Generate all econ summaries, grupped by buidings and fuel prices
tariff_map = uint8([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4; 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]).'; %B, D

RowNames = {'Baseline Costs', 'MGT Costs', 'Savings'};
WriteRange = {'C3', 'C9', 'C15'};
for i = 1:3:36
  f_ind = fuel_index(i);
  d_index = mod(i, 12) + (mod(i, 12) == 0) * 12;
  bd = tariff_map(d_index, 1);
  pause(0.5);
  while (canSave == 0)
    try
      xlswrite(filename, {['Econ Summary, Building ', num2str(bd), ' FC = ', num2str(price_ft3(f_ind)), ' $/1000 ft3']}, ['B', num2str(bd), 'FC', num2str(f_ind)], 'A2');
      canSave = 1;
    catch
    end
  end
  canSave = 0;
  for day = 1:3
    pause(0.5);
    econ_sum = table;
    econ_sum.Electricity = [base_electricity_charge(i+day-1); bought_elec(i+day-1) - sold_energy(i+day-1); 0];
    econ_sum.Heat = [base_heat_charge(i+day-1); bought_heat(i+day-1); 0];
    econ_sum.Fuel = [0; bought_fuel(i+day-1); 0];
    econ_sum.Total = [total_charge(i+day-1); MGT_cost(i+day-1); savings(i+day-1)];
    econ_sum.Properties.RowNames = RowNames;
    while (canSave == 0)
      try
        writetable(econ_sum, filename, 'Sheet', ['B', num2str(bd), 'FC', num2str(f_ind)], 'Range', WriteRange{day}, 'WriteRowNames', true);
        canSave = 1;
      catch
      end
    end
    canSave = 0;
  end
end

%
% Generate all Demand Pricing summaries, grupped by buidings and fuel prices
tariff_map = uint8([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4; 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]).'; %B, D

RowNames = {'PDC', 'IDC'};
WriteRange = {'J3', 'J10', 'J16'};
for i = 1:3:36
  f_ind = fuel_index(i);
  d_index = mod(i, 12) + (mod(i, 12) == 0) * 12;
  bd = tariff_map(d_index, 1);
  pause(0.5);
  while (canSave == 0)
    try
      xlswrite(filename, {['Demand Pricing Summary, Building ', num2str(bd), ' FC = ', num2str(price_ft3(f_ind)), ' $/1000 ft3']}, ['B', num2str(bd), 'FC', num2str(f_ind)], 'J2');
      canSave = 1;
    catch
    end
  end
  canSave = 0;
  for day = 1:3
    demand_sum = table;
    demand_sum.BaseDC = [ut_PDC(i+day-1); ut_IDC(i+day-1)];
    demand_sum.MGTDC = [MGT_PDC(i+day-1); MGT_IDC(i+day-1)];
    demand_sum.Savings = demand_sum.BaseDC - demand_sum.MGTDC;
    demand_sum.Properties.RowNames = RowNames;
    while (canSave == 0)
      try
        writetable(demand_sum, filename, 'Sheet', ['B', num2str(bd), 'FC', num2str(f_ind)], 'Range', WriteRange{day}, 'WriteRowNames', true);
        canSave = 1;
      catch
      end
    end
    canSave = 0;
  end
end
