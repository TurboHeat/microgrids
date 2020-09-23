% TODO: Assign outputs!
function [] = robustShortestPath(timeStepSize, endTime, kwargs)
% Analysis of the solution found by the shortestpath solver
% and comparison to the robust shortest path solvers. 
%
% Shortest path solvers are used according to algorithms in the paper [].
% We look for the path which minimizes the worst-case cost over all
% possible demand profiles modeled. As the price is higher when the
% demand is higher, the "worst-case" considered by the algorithm is
% achieved when the demand is the maximal possible. For that reason,
% we only consider the case in which the demand is equal to 
% "mean + alpha*std", and not the case in which it is equal to
% "mean - alpha*std".
%% Handling inputs:
arguments
  timeStepSize (1,1) double {mustBePositive} = 15 % length of a timestep in [s]
  endTime (1,1) double {mustBePositive} = 24; % final time in [h]
  kwargs.alphaLin (1,1) double = 1; 
  % Parameter for first uncertainty set. 
  % In equation (6) in the paper, we take:
  %   Delta_P(t) = alpha_Linfty*std_P(t) 
  %   Delta_H(t) = alpha_Linfty*std_H(t).
  
  kwargs.alphaMixed (1,1) double = 0.1; 
  % Parameter for second uncertainty set (Algorithm 1 in the paper).
  
  kwargs.alphaSpikes (1,1) double = 1*sqrt( getNumTimesteps(timeStepSize, endTime) ); 
  % Parameter for second uncertainty set (Algorithm 1 in paper).
  % In equation (7) in the paper, we take Delta_P(t) = alpha_Mixed*std_P(t) and
  % Delta_H(t) = alpha_Mixed*std_H(t). Moreover, we take delta_{P,t} = std_P(t),
  % delta_{H,t} = std_H(t) and mu_1 = alpha_spikes.

  kwargs.approx (1,1) ApproximationType = ApproximationType.Exact;
  kwargs.epsilon (1,1) double = 1E-2;
  kwargs.maxIter (1,1) double {mustBeGreaterThanOrEqual(kwargs.maxIter, 2)} = 1E2;
  
  kwargs.dataPath (1,1) string = "../Data"
  kwargs.showTariffs (1,1) logical = false % will plots of tariffs be shown?
  kwargs.showDemands (1,1) logical = false % will plots of demands be shown?
  kwargs.transitionPenalty (1,1) double = 0.01;
end
% Unpack kwargs:
showTariffs = kwargs.showTariffs;
showDemands = kwargs.showDemands;
transitionPenalty = kwargs.transitionPenalty;
epsilon = kwargs.epsilon;

%% Constants:
FUEL_INDEX = repelem(1:3, 1, 12).';
[PRICE_kg_f, HEAT_TARIFF] = NATURAL_GAS_PARAMS();
SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
NO_ENVELOPE_ALPHA = 0;
NODES_CONNECTED_TO_ARTIFICIAL_START = 1;
NO_STATE_TRANSITION_PENALTY = 0;
NO_DEMAND = 0;
% Buildings
BUILDINGS = BuildingType(1:4);
NUM_BUILDINGS = numel(BUILDINGS);

%% Initialization
stepsPerHour = SECONDS_PER_HOUR / timeStepSize; % number of time steps in 1h
%% Loading
tmp = load(fullfile(kwargs.dataPath, 'graph_24h.mat'));
g = tmp.g; SV2State = tmp.svToStateNumber; clear tmp;
[elecTariffs, demands, NUM_WINDOWS] = getTariffsAndDemands(BUILDINGS);
nTotalNodes = size(SV2State, 1); %(smax-smin)*(vmax-vmin+1)+1
% Engine maps:
[POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU] = loadEngineMaps(nTotalNodes, kwargs.dataPath);
% Visualizations?
if showTariffs || showDemands
  [dailyTariffs,tariffQueryTimes] = getDailyTariffs(elecTariffs, timeStepSize); %[nTimestepsPerDay, Elec (1), nDays, nBuildings]
  if showTariffs
    visualizeTariffs(dailyTariffs, stepsPerHour);  
  end
  if showDemands
    upsampledDemands = upsampleDemands(demands, tariffQueryTimes, smoothPeriod); %[nTimestepsPerDay, Elec+Heat (2), nDays, nBuildings, Mean+Std (2)]
    visualizeDemands(upsampledDemands, stepsPerHour);
  end
end

%% Preliminary computations
T = getNumTimesteps(timeStepSize, endTime);
nF = numel(FUEL_INDEX);

[SV_states, timeFrom, nTimesteps, stateFromMap, stateToMap] = ...
  transformAdjacencyMatrix(g.Edges.EndNodes(:,1), g.Edges.EndNodes(:,2), SV2State);

% Select the correct column with the final price
sol_select = [~SV_states(stateFromMap, 1) & ~SV_states(stateToMap, 1), ... % Off-off
  ~SV_states(stateFromMap, 1) & SV_states(stateToMap, 1), ... % Start up
  ~SV_states(stateToMap, 1), ... % Shut down
  true(numel(nTimesteps), 1)]; % Remaining transitions
[~, sol_select] = max(sol_select, [], 2);

% Assign a small penalty to every input (s,v) change
transitionPenaltyFlag = [zeros(nTotalNodes, 1); ...
  ~(SV_states(stateFromMap(nTotalNodes+1:end-nTotalNodes), 1) == SV_states(stateToMap(nTotalNodes+1:end-nTotalNodes), 1) & ... %checks equality of S values
    SV_states(stateFromMap(nTotalNodes+1:end-nTotalNodes), 2) == SV_states(stateToMap(nTotalNodes+1:end-nTotalNodes), 2));
   zeros(nTotalNodes,1)]; % checks equality of V values

% Preallocation
power_nom_MGT = zeros(T, nF);
heat_nom_MGT = zeros(T, nF);
mdot_nom_MGT = zeros(T, nF);
path_nom_MGT = cell(1,nF);
path_cost_nom = zeros(1, nF);
path_edge_nom = cell(1,nF);

power_robust_Linfty_MGT = zeros(T, nF);
heat_robust_Linfty_MGT = zeros(T, nF);
mdot_robust_Linfty_MGT = zeros(T, nF);
path_robust_Linfty = cell(1,nF);
path_cost_robust_Linfty = zeros(1, nF);
path_edge_robust_Linfty = cell(1,nF);

power_robust_Mixed_MGT = zeros(T, nF);
heat_robust_Mixed_MGT = zeros(T, nF);
mdot_robust_Mixed_MGT = zeros(T, nF);
path_robust_Mixed = cell(1,nF);
path_cost_robust_Mixed = zeros(1, nF);
path_edge_robust_Mixed = cell(1,nF);

% Cost of schedule when the true demand is revealed.
nom_cost = zeros(1,nF);
robust_Linfty_cost = zeros(1,nF);
robust_mixed_cost = zeros(1,nF);

%% Compute decided_costs for all algorithms
nPrices = numel(PRICE_kg_f);
progressbar('Gas prices', 'Building types', 'Days');
for iP = 1:nPrices
  fuelPrice = PRICE_kg_f(iP);
  heatTariff = HEAT_TARIFF(iP);
  for iB = 1:NUM_BUILDINGS
    for iW = 1:NUM_WINDOWS      
      d = demands(iW, iB);
      mElec = d.valMean(:,1);
      mHeat = d.valMean(:,2);
      sElec = d.valStd(:,1);
      sHeat = d.valStd(:,2);
            
      decided_costs_nominal = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, NO_ENVELOPE_ALPHA, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);

      decided_costs_robust_Linfty = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alphaLin, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
      
      decided_costs_robust_mixed_nospike = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alphaMixed, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
      
      decided_costs_robust_mixed_withspike = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alphaMixed + kwargs.alphaSpikes, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
      
      %% Find the shortest paths:
      % Run the nominal algorithm
      g_nominal = digraph(stateFromMap, stateToMap, decided_costs_nominal(:, iP));
      [path_nom_MGT{iP}, path_cost_nom(iP), path_edge_nom{iP}] = shortestpath(g_nominal, 1, max(stateToMap), 'Method', 'acyclic');
      [power_nom_MGT(:, iP), heat_nom_MGT(:, iP), mdot_nom_MGT(:, iP)] = extractPath(path_nom_MGT{iP}, power_map, heat_map, fuel_map, SV_states);
      % Run the robust algorithm with Linfty uncertainty set
      g_robust_Linfty = digraph(stateFromMap, stateToMap, decided_costs_robust_Linfty(:, iP));
      [path_robust_Linfty{iP}, path_cost_robust_Linfty(iP), path_edge_robust_Linfty{iP}] = shortestpath(g_robust_Linfty, 1, max(stateToMap), 'Method', 'acyclic');
      [power_robust_Linfty_MGT(:, iP), heat_robust_Linfty_MGT(:, iP), mdot_robust_Linfty_MGT(:, iP)] = extractPath(path_robust_Linfty{iP}, power_map, heat_map, fuel_map, SV_states);
      % Run the robust algorithm with mixed uncertainty set
      W_spike = decided_costs_robust_mixed_withspike(:,iP) - decided_costs_robust_mixed_nospike(:,iP);
      max_W_spike = max(W_spike);
      min_W_spike = min(W_spike);
      unique_W_spike = unique(W_spike);
      nUWS = numel(unique_W_spike);

      if (kwargs.approx == ApproximationType.MaxIter)
          epsilon = (max_W_spike - min_W_spike) / (kwargs.maxIter - 1);
      end

      switch kwargs.approx
        case ApproximationType.Exact
          Thresholds = unique_W_spike;
        case {ApproximationType.Additive, ApproximationType.MaxIter}
          Thresholds = [min_W_spike:epsilon:max_W_spike,max_W_spike].'; %Add max(W_spike) at the end to assure that the maximum is also taken into account
        case ApproximationType.Multiplicative
          Thresholds = exp([log(min_W_spike):log(1+epsilon):log(max_W_spike),log(max_W_spike)]).'; %Add max(W_spike) as before. This is a geometric series.
        otherwise
          error('Unsupported Approximation Type');
      end
      nTrs = numel(Thresholds);
      if nTrs > nUWS
          Thresholds = unique_W_spike;
      end 
      V_costs = zeros(nTrs,1);
      V_paths = cell(nTrs,1);
      V_edge_paths = cell(nTrs,1);
      g_Mixed = digraph(stateFromMap, stateToMap, decided_costs_robust_mixed_nospike(:, iP));
      % Retrieve the index list of the edges - saves time later when changing weights.
      edgePermutationMap = g.Edges.Weight;
      for iT = 1:nTrs
          Weights = (W_spike <= Thresholds(iT)).*decided_costs_robust_mixed_nospike(:, iP) + (W_spike > Thresholds(iT)).*inf;
          g_Mixed.Edges.Weight = Weights(edgePermutationMap);
          [V_paths{iT}, path_cost, edge_path] = shortestpath(g_Mixed, 1, max(stateToMap), 'Method', 'acyclic');
          V_costs(iT) = path_cost + max(W_spike(edge_path));
          V_edge_paths{iT} = edge_path;
      end

      [path_cost_robust_Mixed(:), index_path_Mixed] = min(V_costs);
      path_robust_Mixed(iP) = V_paths(index_path_Mixed);
      path_edge_robust_Mixed{iP} = V_edge_paths(index_path_Mixed);
      [power_robust_Mixed_MGT(:, iP), heat_robust_Mixed_MGT(:, iP), mdot_robust_Mixed_MGT(:, iP)] = extractPath(path_robust_Mixed{iP}, POWER_MAP, HEAT_MAP, FUEL_MAP, SV_states);
      
      %% Compare schedules for a single scenario
      % Replace old cost calculation method, which uses a ZOH for the demands
      % and generations, with a new edge-based method that linearly
      % interpolates the demand and generation.

      %%%TODO: If needed (ask Beni) - separate costs to bought electricity,
      %%% sold electricity, bought fuel, and bought heat. Do this by looking at
      %%% the individual components which generate the edge cost in
      %%% assignCosts.m
      
      pen = path_edge_nom{iP};
      nom_cost(iP) = sum(assignCostsInternal(...
        sol_select(pen), stateFromMap(pen), stateToMap(pen), NODES_CONNECTED_TO_ARTIFICIAL_START, ...
        NO_DEMAND, NO_DEMAND, NO_DEMAND, NO_DEMAND, NO_ENVELOPE_ALPHA, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag(pen), NO_STATE_TRANSITION_PENALTY, ...
        timeFrom(pen), nTimesteps(pen), stepsPerHour, timeStepSize));
      
      pel = path_edge_robust_Linfty{iP};
      robust_Linfty_cost(iP) = sum(assignCostsInternal(...
        sol_select(pel), stateFromMap(pel), stateToMap(pel), NODES_CONNECTED_TO_ARTIFICIAL_START, ...
        NO_DEMAND, NO_DEMAND, NO_DEMAND, NO_DEMAND, kwargs.alphaLin, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag(pel), NO_STATE_TRANSITION_PENALTY, ...
        timeFrom(pel), nTimesteps(pel), stepsPerHour, timeStepSize));      
      
      pem = path_edge_robust_Mixed{iP};
      robust_mixed_cost(iP) = sum(assignCostsInternal(...
        sol_select(pem), stateFromMap(pem), stateToMap(pem), NODES_CONNECTED_TO_ARTIFICIAL_START, ...
        NO_DEMAND, NO_DEMAND, NO_DEMAND, NO_DEMAND, kwargs.alphaMixed, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag(pem), NO_STATE_TRANSITION_PENALTY, ...
        timeFrom(pem), nTimesteps(pem), stepsPerHour, timeStepSize));
      
      % TODO: Make sure the above are equivalent to:      
      %{
      nom_cost(iP) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                        sol_select(path_edge_nom), time_from(path_edge_nom), n_tsteps(path_edge_nom),...
                        from_state_map(path_edge_nom), to_state_map(path_edge_nom),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                        elec_tariff(:, d_index), heat_tariff(fuel_index(iP)), price_kg_f(fuel_index(iP)), transition_penalty(path_edge_nom),0));  

      robust_Linfty_cost(iP) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                        sol_select(path_edge_robust_Linfty), time_from(path_edge_robust_Linfty), n_tsteps(path_edge_robust_Linfty),...
                        from_state_map(path_edge_robust_Linfty), to_state_map(path_edge_robust_Linfty),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                        elec_tariff(:, d_index), heat_tariff(fuel_index(iP)), price_kg_f(fuel_index(iP)), transition_penalty(path_edge_robust_Linfty),0));  

      robust_mixed_cost(iP) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                        sol_select(path_edge_robust_Mixed), time_from(path_edge_robust_Mixed), n_tsteps(path_edge_robust_Mixed),...
                        from_state_map(path_edge_robust_Mixed), to_state_map(path_edge_robust_Mixed),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                        elec_tariff(:, d_index), heat_tariff(fuel_index(iP)), price_kg_f(fuel_index(iP)), transition_penalty(path_edge_robust_Mixed),0));
      %}              
      keyboard; % TODO: make sure some aggregation (iteration results are saved somewhere) is taking place!
    end
  end
end

end

function [decided_costs] = assignCostsInternal(...
  sol_select, stateFromMap, stateToMap, nTotalNodes, ...
  elecDemandMean, elecDemandStd, heatDemandMean, heatDemandStd, demandStandardEnvelope, ...
  elecTariff, heatTariff, fuelPrice,...
  powerMap, heatMap, fuelMap, mdot_fuel_SU, ...
  transitionPenaltyFlag, transitionPenalty, ...
  time_from, nTimesteps, stepsPerHour, timeStepSize)

% Apply alpha:
elecDemand = elecDemandMean + demandStandardEnvelope * elecDemandStd;
heatDemand = heatDemandMean + demandStandardEnvelope * heatDemandStd;

% Assign edge costs
decided_costs = assignCosts(...
  timeStepSize, powerMap, heatMap, fuelMap, ...
  mdot_fuel_SU, nTotalNodes, sol_select, time_from, nTimesteps,...
  stateFromMap, stateToMap,...
  elecDemand, heatDemand, elecTariff, heatTariff, fuelPrice, ...
  transitionPenaltyFlag, transitionPenalty, stepsPerHour);
end

function T = getNumTimesteps(dt, endTime)
  SECONDS_PER_MINUTE = 60;
  MINUTES_PER_HOUR = 60;
  SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
  T = endTime * SECONDS_PER_HOUR / dt;    
end

function [elecTariffs, demands, NUM_WINDOWS] = getTariffsAndDemands(buildings)
[CHP, NUM_WINDOWS] = LOAD_DEMAND_DATASETS();

%% Get all tariff & demand combinations
nBuildings = numel(buildings);
elecTariffs = cell(nBuildings, NUM_WINDOWS);
demands = cell(nBuildings, NUM_WINDOWS);

parfor b = 1:nBuildings
  chp = CHP(b);
  for d = 1:NUM_WINDOWS
    do = chp.next(); demands{b,d} = do; % demand object
    elecTariffs{b,d} = getSeasonalElectricityTariff(b, do.timeEnd);
  end
end
% Unbox:
demands = reshape([demands{:}], nBuildings, NUM_WINDOWS).';
elecTariffs = reshape([elecTariffs{:}], nBuildings, NUM_WINDOWS).';
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

function [powerMap, heatMap, fuelMap, mdotFuelSU] = loadEngineMaps(expectedTotalNodes, dataPath)
arguments 
  expectedTotalNodes (1,1) double {mustBeInteger, mustBePositive} = NaN
  dataPath (1,1) string = "../Data"
end
%% Load (and optionally rename) mappings
if isfile(fullfile(dataPath, 'CHP.mat')) % Engine similar to Capstone C65 (incl. extrapolation)
  tmp = load(fullfile(dataPath, 'CHP.mat'));
  powerMap = tmp.PowerSurf;
  heatMap = tmp.PheatSurfExt;
  m_dot_f_map = tmp.FFSurfExt;
else
  tmp = load(fullfile(dataPath, 'sv_mappings.mat')); % "Virtual" engine
  powerMap = tmp.power_map;
  heatMap = tmp.heat_map;
  m_dot_f_map = tmp.m_dot_f_map;
end
% Verify that we have the right mappings:
assert( isequal( numel(powerMap), numel(heatMap), numel(m_dot_f_map), expectedTotalNodes-1 ), ...
  'Mismatch between matrix sizes! Try using a different ''sv_mappings'' file.');

% Reshape the power maps
powerMap = [reshape(powerMap.', [], 1); 0];
heatMap = [reshape(heatMap.', [], 1); 0];
fuelMap = [reshape(m_dot_f_map.', [], 1); 0];
mdotFuelSU = sum(m_dot_f_map(:, 1:end-1)+m_dot_f_map(:, 2:end), 2); % Computes mdot of fuel consumed in startup for all the different max states
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
  [nTimestepsPerDay, ~, nDays, nBuildings] = size(t); %#ok<ASGLU>
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