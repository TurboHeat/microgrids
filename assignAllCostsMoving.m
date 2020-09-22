function [varargout] = assignAllCostsMoving(kwargs)
%% Handle inputs
arguments
  kwargs.showTariffs (1,1) logical = false % will plots of tariffs be shown?
  kwargs.showDemands (1,1) logical = false % will plots of demands be shown?
  kwargs.smoothDemandTimesteps (1,1) double {mustBeInteger, mustBePositive} = 21; % number of timesteps for smoothing. 1=off
  kwargs.endTime (1,1) double {mustBePositive} = 24; % [h]
  kwargs.savePath (1,1) string = "../Data/";
  kwargs.transitionPenalty (1,1) double = 0.01;
  kwargs.demandStandardEnvelope (1,1) double {mustBeNonnegative} = 0; % α in: expectedDemand = μ ± α·σ
  kwargs.timeStepSize (1,1) double = 15; % duration of time step in [s]
end
% Unpack kwargs:
showTariffs = kwargs.showTariffs;
showDemands = kwargs.showDemands;
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
NUM_BUILDINGS = numel(BUILDING);

% Demands
[CHP, NUM_WINDOWS] = LOAD_DEMAND_DATASETS();

%% Get all tariff & demand combinations
elecTariffs = cell(NUM_BUILDINGS, NUM_WINDOWS);
demands = cell(NUM_BUILDINGS, NUM_WINDOWS);

parfor b = 1:NUM_BUILDINGS
  chp = CHP(b);
  for d = 1:NUM_WINDOWS
    do = chp.next(); demands{b,d} = do; % demand object
    elecTariffs{b,d} = getSeasonalElectricityTariff(b, do.timeEnd);
  end
end
% Unbox:
demands = reshape([demands{:}], NUM_BUILDINGS, NUM_WINDOWS).';
elecTariffs = reshape([elecTariffs{:}], NUM_BUILDINGS, NUM_WINDOWS).';

%% Visualizations?
if showTariffs || showDemands
  [dailyTariffs,tariffQueryTimes] = getDailyTariffs(elecTariffs, dt); %[nTimestepsPerDay, Elec (1), nDays, nBuildings]
end
if showTariffs
  % Get tariff at each time step:
  visualizeTariffs(dailyTariffs, N_LINES);  
end
if showDemands
  upsampledDemands = upsampleDemands(demands, tariffQueryTimes, smoothPeriod); %[nTimestepsPerDay, Elec+Heat (2), nDays, nBuildings, Mean+Std (2)]
  visualizeDemands(upsampledDemands, N_LINES);
end

%% Run Iliya's mapping and save all the variables at once
load(fullfile(savePath, 'graph_24h.mat'), 'g', 'svToStateNumber');
state_from = g.Edges.EndNodes(:,1);
state_to = g.Edges.EndNodes(:,2);

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
nPrices = numel(PRICE_kg_f);
% decided_costs = zeros(nStates, nPrices, NUM_BUILDINGS, NUM_WINDOWS);
progressbar('Gas prices', 'Building types', 'Days');
for iP = 1:nPrices
  fuel_price = PRICE_kg_f(iP);
  for iB = 1:NUM_BUILDINGS
    for iW = 1:NUM_WINDOWS
      d = demands(iW, iB);
      mPower = d.valMean(:,1);
      mHeat = d.valMean(:,2);
      sPower = d.valStd(:,1);
      sHeat = d.valStd(:,2);
      decided_costs = assignCosts(...
        double(dt), power_map, heat_map, fuel_map, ...
        mdot_fuel_SU, total_nodes, sol_select, time_from, n_tsteps,...
        from_state_map, to_state_map,...
        mPower + alpha * sPower, mHeat + alpha * sHeat, ...
        elecTariffs(iW, iB), HEAT_TARIFF(iP), fuel_price, ...
        transition_penalty_indicator, transitionPenalty, N_LINES);
      % TODO: do something with decided_costs, or save a summary
      progressbar([], [], iW/NUM_WINDOWS);
      keyboard; % decided_costs isn't being processed in any way!!!
    end
    progressbar([], iB/NUM_BUILDINGS, []);
  end
  progressbar(iP/nPrices, [], []);
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
% Count how many times "next" can be called:
nextCnt = 0; % Should be equal to: 365-14+1 = 352
while (chp(1).hasNext)
  nextCnt = nextCnt + 1;
  [~] = chp(1).next();
end
% Rewind again:
arrayfun(@(x)x.fastForward(startDay, 'last'), chp);
%}
chp = struct2array(load('../Data/CHP2004.mat', 'chp'));
nWindows = struct2array(load('../Data/CHP2004.mat', 'nextCnt'));
end

function [dailyTariffs,tariffQueryTimes] = getDailyTariffs(hourlyElecTariffs, dt)
SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
HOURS_PER_DAY = 24;
stepsPerHour = SECONDS_PER_HOUR / dt; %number of lines per hour
tariffQueryTimes = linspace(0, HOURS_PER_DAY, stepsPerHour*HOURS_PER_DAY+1).'; tariffQueryTimes(end) = [];
dailyTariffs = cell2mat(shiftdim(...
  arrayfun( @(x)tariffAtTime(x, tariffQueryTimes), hourlyElecTariffs, 'UniformOutput', false), -2));
end

function usDemand = upsampleDemands(demands, queryToD, smoothingWindowSz)
arguments
  demands (:,:) TimestampedStatHolder
  queryToD (:,1) double
  smoothingWindowSz (1,1) double {mustBeInteger, mustBeNonnegative}
end
% This function returns a 5D array where:
% Dim1 is the time-of-day
% Dim2 is either electricity (at index 1) or heat (at index 2)
% Dim3 is the day-of-year
% Dim4 is the building type
% Dim5 is either mean (at index 1) or stdev (at index 2)
%% Constants:
COLUMNS_FOR_POWER_DATA = 1;
COLUMNS_FOR_HEAT_DATA = 1;
nC = COLUMNS_FOR_POWER_DATA + COLUMNS_FOR_HEAT_DATA;
%% Input Handling:
X = demands(1).hourOfVal;
nT = numel(queryToD);
[nW, nB] = size(demands); % NUM_WINDOWS, NUM_BUILDINGS
% Preallocation:
[meanDemand,stdevDemand] = deal(zeros(nT, nC, nW, nB));

% Interpolation of the means and standard deviations:
tmpM = cell2mat(reshape({demands.valMean},1,1,nW,nB));
tmpS = cell2mat(reshape({demands.valStd},1,1,nW,nB));
for indQ = 1:nC
  for indD = 1:nW
    parfor indB = 1:nB
      meanDemand(:, indQ, indD, indB)  = nakeinterp1(X, tmpM(:,indQ, indD, indB), queryToD);
      stdevDemand(:, indQ, indD, indB) = nakeinterp1(X, tmpS(:,indQ, indD, indB), queryToD);
    end
  end
end

% Time smoothing (using convolution):
if smoothingWindowSz > 1
  convKer = ones(smoothingWindowSz, 1, 1, 1, 1); % The only non-singleton dimension corresponds to time.
  usDemand = convn(cat(5, meanDemand, stdevDemand), convKer, 'same');
else
  usDemand = cat(5, meanDemand, stdevDemand);
end
end

function [] = visualizeTariffs(t, stepsPerHour)
  [nTimestepsPerDay, ~, nDays, nBuildings] = size(t);
  assert(nBuildings == 4, 'The function is configured to visualize exactly 4 buildings.');
  %% Preparations    
  BLDG_DESCR = {'Large Hotel', 'Full Service Restaraunt', 'Small Hotel', 'Residential'};
  tsId = 1:nTimestepsPerDay;
  hF = figure();
  hTL = tiledlayout(hF, 2,2,'TileSpacing','compact');
  hAx = gobjects(2,2);
  t1 = permute(t, [3,1,4,2]);
  %% 3D Version:
  %{
  DAYS_IN_TWO_WEEKS = 14;
  dayId = (1:nDays)+DAYS_IN_TWO_WEEKS;
  [TT,DD] = meshgrid(tsId/stepsPerHour, dayId);
  for iB = 1:nBuildings %building
    hAx(iB) = nexttile(hTL, iB);
    waterfall(hAx(iB), DD, TT, t1(:,:,iB));
    title(hAx(iB), BLDG_DESCR{iB});
  end
  set(hAx, 'YDir', 'reverse', 'YLim', [0 24], 'YTick', 0:4:24, ...
           'XLim', [14 366], 'XTick', 14:44:366);
  xlabel(hAx, 'Day of Year')
  ylabel(hAx, 'Time of Day [h]');
  zlabel(hAx, 'Rate [$/kWh]')  
  %}
  %% 2D Version:
  NON_SUMMER_DAY_ID = 1;
  SUMMER_DAY_ID = 200;
  for iB = 1:nBuildings %building
    hAx(iB) = nexttile(hTL, iB);
    reset(hAx(iB));
    hL = plot(hAx(iB), ...
      tsId/stepsPerHour, t1(NON_SUMMER_DAY_ID,:,iB),...
      tsId/stepsPerHour, t1(SUMMER_DAY_ID,:,iB), '--', 'LineWidth', 2);
    hL(1).DisplayName = '"Winter" tariff';
    hL(2).DisplayName = '"Summer" tariff';
    title(hAx(iB), BLDG_DESCR{iB});
    legend(hAx(iB), 'Location', 'northwest');
  end
  set(hAx, 'XLim', [0 24], 'XTick', 0:2:24, 'FontSize', 16); 
  grid(hAx,'on');
  xlabel(hAx, 'Time of Day [h]');
  ylabel(hAx, 'Rate [$/kWh]');
end

function [] = visualizeDemands(d, stepsPerHour)
  [nTimestepsPerDay, ~, nDays, nBuildings,~] = size(d); %[nTimestepsPerDay, Elec+Heat (2), nDays, nBuildings, Mean+Std (2)]
  assert(nBuildings == 4, 'The function is configured to visualize exactly 4 buildings.');
  %% Preparations    
  BLDG_DESCR = {'Large Hotel', 'Full Service Restaraunt', 'Small Hotel', 'Residential'};
  tsId = 1:nTimestepsPerDay;
  hF = figure();
  hTL = tiledlayout(hF, 2, 4, 'TileSpacing', 'compact');
  hAx = gobjects(2,4);
  md = permute(d(:,:,:,:,1), [3,1,4,2]);
  %% 3D Version:
  DAYS_IN_TWO_WEEKS = 14;
  dayId = (1:nDays)+DAYS_IN_TWO_WEEKS;
  [TT,DD] = meshgrid(tsId/stepsPerHour, dayId);
  for iB = 1:nBuildings %building
    % Electricity
    hAx(iB) = nexttile(hTL, iB);
    waterfall(hAx(iB), DD, TT, md(:,:,iB,1));
    title(hAx(iB), BLDG_DESCR{iB});
    % Power
    hAx(iB+4) = nexttile(hTL, iB+4);
    waterfall(hAx(iB+4), DD, TT, md(:,:,iB,2));
  end
  set(hAx, 'YDir', 'reverse', 'YLim', [0 24], 'YTick', 0:4:24, ...
           'XLim', [14 366], 'XTick', 14:44:366);
  xlabel(hAx, 'Day of Year')
  ylabel(hAx, 'Time of Day [h]');
  zlabel(hAx(1:4), 'Electric Demand [kWh]')
  zlabel(hAx(5:8), 'Heat Demand [kWh]')
end