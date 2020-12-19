function [statistics] = analyzeResults(resultsFolder)
arguments
  resultsFolder (1,1) string = '../Data/Results/'
end

mergedFN = fullfile(resultsFolder, "all_results.mat");
if isfile(mergedFN)
  cases = struct2array( load(mergedFN) );
else
%% Build a metadata matrix:
% First use the filenames
[cases] = parseResultFilenames(resultsFolder);
nCases = size(cases,1);
% Reformat some fields:
cases.("Building Type") = BuildingType(cases.("Building Type"));
% Add some derived variables:
cases.("Was Benchmark Run?") = cases.("Algorithm Type") == -1;
cases.("Was Nominal Run?") = cases.("Algorithm Type") == 0;
cases.("Was Robust Run?") = cases.("Algorithm Type") == 1;
cases.("Was Mixed Run?") = cases.("Algorithm Type") == 2;
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
save(mergedFN, 'cases');

%% Cleanup:
clearvars('-except', 'resultsFolder', 'cases') 
end

%% TODO: Integrate
%{
% Initialize Structure - the values have no meaning, only the data
% structure.
currentScenarioData(numAlgorithms) = outputData(1).Data;
allScenariosData(numScenarios).Data = currentScenarioData;
allScenariosData(numScenarios).BuildingType = 0;
allScenariosData(numScenarios).PriceIndex = 0;

for jScenario = 1:numScenarios
  building = buildingTypes(jScenario);
  priceInd = priceIndices(jScenario);
  k = 1;
  for iter = 1:numPairs
    if (building == outputData(iter).BuildingType) && ...
       (priceInd == outputData(iter).PriceIndex)
      currentScenarioData(k) = outputData(iter).Data;
      k = k+1;
    end
  end
  allScenariosData(jScenario).Data = currentScenarioData;
  allScenariosData(jScenario).BuildingType = building;
  allScenariosData(jScenario).PriceIndex = priceInd;
end
%}

% Build a structure to store statistics
statistics.PerDay(nCases).TrueCost = [];
statistics.PerDay(nCases).TrueOverEstimateCost = [];
if (1 == WasBenchmarkRun)
  statistics.PerDay(nCases).CostOverOptimal = [];
end
if (1 == WasNominalRun)
  statistics.PerDay(nCases).DifferenceFromNominal = [];
end

statistics.PerYear(nCases).TrueCost = [];
statistics.PerYear(nCases).TrueOverEstimateCost = [];
if (1 == WasBenchmarkRun)
  statistics.PerYear(nCases).CostOverOptimal = [];
end
if (1 == WasNominalRun)
  statistics.PerYear(nCases).DifferenceFromNominal = [];
end

for id = 1:nCases
  CurrentScenarioData = AllScenariosData(id).Data;
  TitleSuffix = [' - Scenario #', num2str(id)];
  
  
  %Plot "true" cost of schedule - cost of schedule for ground-truth
  %demand
  figure(100*(id - 1)+1);
  for j = 1:length(CurrentScenarioData)
    plot(CurrentScenarioData(j).TrueCost);
    hold all;
  end
  grid on;
  xlabel('Day in year (from January 15^{th}');
  ylabel('Cost [$]');
  title(['Cost for Ground Truth Demand', TitleSuffix]);
  legend(Legend);
  
  
  %Plot Ratio of Estimated cost over True Cost (Budget exceeding)
  figure(100*(id - 1)+2);
  for j = 2:length(CurrentScenarioData)
    plot((CurrentScenarioData(j).TrueCost ./ CurrentScenarioData(j).EstimatedCost - 1)*100);
    hold all;
  end
  grid on;
  xlabel('Day in year (from January 15^{th}');
  ylabel('Cost over Expected [%]');
  title(['Cost over Estimate ("Over Budget")', TitleSuffix]);
  legend(Legend(2:end));
  
  for j = 1:length(CurrentScenarioData)
    TrueCost = CurrentScenarioData(j).TrueCost;
    
    statistics.PerDay(id).TrueCost = TrueCost;
    statistics.PerYear(id).TrueCost(j).Average = mean(TrueCost);
    statistics.PerYear(id).TrueCost(j).Std = std(TrueCost);
    statistics.PerYear(id).TrueCost(j).Max = max(TrueCost);
    statistics.PerYear(id).TrueCost(j).Min = min(TrueCost);
    statistics.PerYear(id).TrueCost(j).Median = median(TrueCost);
    
    Ratio = (CurrentScenarioData(j).TrueCost ./ CurrentScenarioData(j).EstimatedCost - 1) * 100;
    statistics.PerDay(id).TrueOverEstimateCost = Ratio;
    statistics.PerYear(id).TrueOverEstimateCost(j).TimePositive = sum(Ratio > 0) / length(Ratio) * 100;
    statistics.PerYear(id).TrueOverEstimateCost(j).Max = max(Ratio);
    statistics.PerYear(id).TrueOverEstimateCost(j).ConditionalValueAtRisk = mean(Ratio(Ratio > 0));
    statistics.PerYear(id).TrueOverEstimateCost(j).StdofValueAtRisk = std(Ratio(Ratio > 0));
  end
  
  if (1 == WasBenchmarkRun)
    BenchmarkData = CurrentScenarioData(1);
    
    %Compute True cost over Nominal True Cost.
    figure(100*(id - 1)+3);
    for j = 2:length(CurrentScenarioData)
      plot((CurrentScenarioData(j).TrueCost ./ BenchmarkData.TrueCost - 1)*100);
      hold all;
    end
    grid on;
    xlabel('Day in year (from January 15^{th}');
    ylabel('Cost over Benchmark [%]');
    title(['Cost over Benchmark for Ground Truth Demand', TitleSuffix]);
    legend(Legend(2:end));
    
    for j = 1:length(CurrentScenarioData)
      Ratio = (CurrentScenarioData(j).TrueCost ./ BenchmarkData.TrueCost - 1) * 100;
      statistics.PerDay(id).CostOverOptimal = Ratio;
      if (j == 1) %Benchmark algorithm
        statistics.PerYear(id).CostOverOptimal(j).TimePositive = 0;
        statistics.PerYear(id).CostOverOptimal(j).Max = 0;
        statistics.PerYear(id).CostOverOptimal(j).ConditionalValueAtRisk = 0;
        statistics.PerYear(id).CostOverOptimal(j).StdofValueAtRisk = 0;
      else %Other algorithms
        statistics.PerYear(id).CostOverOptimal(j).TimePositive = sum(Ratio > 0) / length(Ratio) * 100;
        statistics.PerYear(id).CostOverOptimal(j).Max = max(Ratio);
        statistics.PerYear(id).CostOverOptimal(j).ConditionalValueAtRisk = mean(Ratio(Ratio > 0));
        statistics.PerYear(id).CostOverOptimal(j).StdofValueAtRisk = std(Ratio(Ratio > 0));
      end
    end
  end
    
  if (1 == WasNominalRun)
    NominalData = CurrentScenarioData(2);
    
    %Compute how different the solution is different from the nominal solution.
    figure(100*(id - 1)+4);
    for j = 3:length(CurrentScenarioData)
      plot(sum((CurrentScenarioData(j).Power_Generation ~= NominalData.Power_Generation | ...
        CurrentScenarioData(j).Heat_Generation ~= NominalData.Heat_Generation | ...
        CurrentScenarioData(j).Fuel_Consumption ~= NominalData.Fuel_Consumption), 2));
      hold all;
    end
    grid on;
    xlabel('Day in year (from January 15^{th}');
    ylabel('#hours with different schedule from nominal');
    title(['Difference from Nominal Solution', TitleSuffix]);
    legend(Legend(3:end));
    
    for j = 1:length(CurrentScenarioData)
      
      Diff = sum((CurrentScenarioData(j).Power_Generation ~= NominalData.Power_Generation | ...
        CurrentScenarioData(j).Heat_Generation ~= NominalData.Heat_Generation | ...
        CurrentScenarioData(j).Fuel_Consumption ~= NominalData.Fuel_Consumption), 2);
      
      statistics.PerDay(id).DifferenceFromNominal = Diff;
      statistics.PerYear(id).DifferenceFromNominal(j).Max = max(Diff);
      statistics.PerYear(id).DifferenceFromNominal(j).Min = min(Diff);
      statistics.PerYear(id).DifferenceFromNominal(j).Average = mean(Diff);
      statistics.PerYear(id).DifferenceFromNominal(j).Median = median(Diff);
    end
  end
end

end

function tab = loadDataFile(path)
tab = struct2table(struct2array(load(path, 'outputData')), 'AsArray', true);
end