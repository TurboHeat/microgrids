function [statistics] = analyzeResults(resultsFolder)
arguments
  resultsFolder (1,1) string = '../Data/Results/'
end
%% Stage 1: Prepare raw data + metadata
mergedFn = fullfile(resultsFolder, "all_results.mat");
if isfile(mergedFn)
  cases = struct2array( load(mergedFn) );
  nCases = size(cases,1);
else
  %% Build a metadata matrix:
  % First use the filenames
  [cases] = parseResultFilenames(resultsFolder);
  nCases = size(cases,1);
  % Reformat some fields:
  cases.("Algorithm Type") = AlgorithmType(cases.("Algorithm Type"));
  cases.("Building Type") = BuildingType(cases.("Building Type"));
  % Add some derived variables:
  cases.("Was Benchmark Run?") = cases.("Algorithm Type") == AlgorithmType.Benchmark;
  cases.("Was Nominal Run?") = cases.("Algorithm Type") == AlgorithmType.Nominal;
  cases.("Was Robust Run?") = cases.("Algorithm Type") == AlgorithmType.L_inf;
  cases.("Was Mixed Run?") = cases.("Algorithm Type") == AlgorithmType.Mixed;
  cases.("Robust Index") = cumsum(cases.("Was Robust Run?")).*cases.("Was Robust Run?");
  cases.("Mixed Index") = cumsum(cases.("Was Mixed Run?")).*cases.("Was Mixed Run?");
  % Create legend entries:
  cases.Legend = strings(nCases,1);
  cases(cases.("Was Benchmark Run?"), "Legend") = table("Benchmark");
  cases(cases.("Was Nominal Run?"),   "Legend") = table("Nominal");
  id = cases.("Was Robust Run?");
  cases(id, "Legend") = table("Rob. L_\infty #" + cases{id, "Robust Index"});
  id = cases.("Was Mixed Run?");
  cases(id, "Legend") = table("Rob. Mixed #" + cases{id, "Mixed Index"});

  %% Load data files:
  for id = nCases:-1:1
    od(id,:) = loadDataFile(fullfile(resultsFolder, cases{id, "name"}));
  end
  % Remove variables:
  od = removevars(od, 'AlgorithmType');
  % Rename variables:
  od.Properties.VariableNames = ["Power Generation", "Heat Generation", ...
    "Fuel Consumption", "Estimated Cost", "True Cost", "Algorithm Parameters"];
  % Merge tables:
  cases = [cases, od];
  % Reorder variables:
  cases = movevars(cases, 'Algorithm Parameters', 'After', "Algorithm Parameters ID");

  %% Save
  save(mergedFn, 'cases');

  %% Cleanup:
  clearvars('-except', 'resultsFolder', 'cases', 'nCases') 
end
% Cleanup
clear mergedFn

%% Stage 2: Derived Statistics
statsFn = fullfile(resultsFolder, 'all_statistics.mat');
if isfile(statsFn)
  tmp = load(statsFn);
  G = tmp.G;
  nG = tmp.nG;
  statistics = tmp.statistics;
  clear tmp
else
  %% Group runs into "scenarios" (unique combinations of building and fuel prices)
  G = findgroups(cases.("Building Type"), cases.("Fuel Price ID"), cases.("Power Scaling Factor"));
  nG = max(G);

  %% Compute statistics:
  % Preallocate a structure to store statistics
  statistics.PerDay(nCases).TrueCost = [];
  statistics.PerDay(nCases).TrueOverEstimateCost = [];
  statistics.PerDay(nCases).CostOverOptimal = [];
  statistics.PerDay(nCases).DifferenceFromNominal = [];
  statistics.PerYear(nCases).TrueCost = [];
  statistics.PerYear(nCases).TrueOverEstimateCost = [];
  statistics.PerYear(nCases).CostOverOptimal = [];
  statistics.PerYear(nCases).DifferenceFromNominal = [];
  % Gather statistics for each scenario
  for indG = 1:nG
    idx = G==indG;
    [statistics.PerDay(idx), statistics.PerYear(idx)] = getStatistics(cases(idx,:));
  end
  save( statsFn, 'G', 'nG', 'statistics' );
end
% Cleanup
clear statsFn

%% Stage 3: Visualization
for indG = 1:nG
  idx = G==indG;
  hF = visualizeScenario(cases(idx,:));
  set(hF,'Name', "Scenario #" + indG);
end

end

function [dayStats, yearStats] = getStatistics(scenario)
%% Other preparations
nRuns = size(scenario, 1);
% Preallocation:
[dayStats, yearStats] = deal(repmat(struct('TrueCost',[],'TrueOverEstimateCost',[],...
  'CostOverOptimal',[],'DifferenceFromNominal',[]),[nRuns,1]));

bc = scenario( scenario.("Was Benchmark Run?"),:);
nc = scenario( scenario.("Was Nominal Run?"),:);
for iR = 1:nRuns
  tc = scenario.("True Cost"){iR};
  ec = scenario.("Estimated Cost"){iR};
  ratio = (tc ./ ec - 1)*100;
  
  dayStats(iR) .TrueCost = tc;
  yearStats(iR).TrueCost.Average = mean(tc);
  yearStats(iR).TrueCost.Std = std(tc);
  yearStats(iR).TrueCost.Max = max(tc);
  yearStats(iR).TrueCost.Min = min(tc);
  yearStats(iR).TrueCost.Median = median(tc);
  
  dayStats(iR) .TrueOverEstimateCost = ratio;
  yearStats(iR).TrueOverEstimateCost.TimePositive = sum(ratio > 0) / length(ratio) * 100;  
  yearStats(iR).TrueOverEstimateCost.Max = max(ratio);
  yearStats(iR).TrueOverEstimateCost.ConditionalValueAtRisk = nanmean(ratio(ratio > 0));
  yearStats(iR).TrueOverEstimateCost.StdofValueAtRisk = nanstd(ratio(ratio > 0));
  
  ratio = (tc ./ bc.("True Cost"){1} - 1) * 100;
  dayStats(iR).CostOverOptimal = ratio;
  if (iR == 1) %Benchmark algorithm
    yearStats(iR).CostOverOptimal.TimePositive = 0;
    yearStats(iR).CostOverOptimal.Max = 0;
    yearStats(iR).CostOverOptimal.ConditionalValueAtRisk = 0;
    yearStats(iR).CostOverOptimal.StdofValueAtRisk = 0;
  else %Other algorithms
    yearStats(iR).CostOverOptimal.TimePositive = sum(ratio > 0) / length(ratio) * 100;
    yearStats(iR).CostOverOptimal.Max = max(ratio);
    yearStats(iR).CostOverOptimal.ConditionalValueAtRisk = mean(ratio(ratio > 0));
    yearStats(iR).CostOverOptimal.StdofValueAtRisk = std(ratio(ratio > 0));
  end  
  
  Diff = sum((scenario.("Power Generation"){iR} ~= nc.("Power Generation"){1} | ...
              scenario.("Heat Generation"){iR}  ~= nc.("Heat Generation"){1} | ...
              scenario.("Fuel Consumption"){iR} ~= nc.("Fuel Consumption"){1}), 2);
  
  dayStats(iR).DifferenceFromNominal = Diff;
  yearStats(iR).DifferenceFromNominal.Max = max(Diff);
  yearStats(iR).DifferenceFromNominal.Min = min(Diff);
  yearStats(iR).DifferenceFromNominal.Average = mean(Diff);
  yearStats(iR).DifferenceFromNominal.Median = median(Diff);  
end

end

function hF = visualizeScenario(scenario)
NUM_FIGURES = 4;
nRuns = size(scenario, 1);
hF = gobjects(NUM_FIGURES,1); hAx = gobjects(NUM_FIGURES,1); 
hL = gobjects(NUM_FIGURES,1);
bc = scenario( scenario.("Was Benchmark Run?"),:);
nc = scenario( scenario.("Was Nominal Run?"),:);
% Create figures:
for k = 1:NUM_FIGURES
  hF(k) = figure('NumberTitle','off'); 
  hAx(k) = axes(hF(k));  %#ok<LAXES>
  hL(k) = legend(hAx(k), 'Location', 'eastoutside', 'NumColumns', 3);  
end
set(hL, 'LimitMaxLegendEntries', false);
hold(hAx, 'on');
grid(hAx,'on');
%% Figure 1: "true" cost of schedule (cost of schedule for ground-truth demand)
for iR = 1:nRuns
  tc = scenario.("True Cost"){iR};
  plot(hAx(1), tc, 'DisplayName', scenario.("Legend")(iR));
end
xlabel(hAx(1),'Day in year (from January 15^{th}');
ylabel(hAx(1),'Cost [$]');
title(hAx(1),'Cost for Ground Truth Demand');

%% Figure 2: ratio of the true cost and the estimated cost ("budget exceeding")
for iR = 1:nRuns
  tc = scenario.("True Cost"){iR};
  ec = scenario.("Estimated Cost"){iR};
  ratio = (ec ./ tc - 1)*100;  
  plot(hAx(2), ratio, 'DisplayName', scenario.("Legend")(iR));
end
xlabel(hAx(2),'Day in year (from January 15^{th}');
ylabel(hAx(2),'Cost over Expected [%]');
title(hAx(2),'Cost over Estimate ("Over Budget")');

%% Figure 3: True cost over the nominal (un-robustified) cost
for iR = 2:nRuns
  tc = scenario.("True Cost"){iR};
  plot(hAx(3), (tc ./ bc.("True Cost"){1} - 1)*100, ...
       'DisplayName', scenario.("Legend")(iR));
end
xlabel(hAx(3), 'Day in year (from January 15^{th}');
ylabel(hAx(3), 'Cost over Benchmark [%]');
title(hAx(3), 'Cost over Benchmark for Ground Truth Demand');

%% Figure 4: Different between solution N and the nominal solution
for iR = 3:nRuns
  Diff = sum((scenario.("Power Generation"){iR} ~= nc.("Power Generation"){1} | ...
              scenario.("Heat Generation"){iR}  ~= nc.("Heat Generation"){1} | ...
              scenario.("Fuel Consumption"){iR} ~= nc.("Fuel Consumption"){1}), 2);
  plot(hAx(4), Diff, 'DisplayName', scenario.("Legend")(iR));
end
xlabel(hAx(4), 'Day in year (from January 15^{th}');
ylabel(hAx(4), '#hours with different schedule from nominal');
title(hAx(4), 'Difference from Nominal Solution');
end

function tab = loadDataFile(path)
tab = struct2table(struct2array(load(path, 'outputData')), 'AsArray', true);
end