function [Statistics] = parseData(AllScenariosData)
close all;


NumOfScenarios = length(AllScenariosData);
NumOfAlgorithms = length(AllScenariosData(1).Data);
WasBenchmarkRun = AllScenariosData(1).Data(1).AlgorithmType == -1;
WasNominalRun = AllScenariosData(1).Data(2).AlgorithmType == 0;


% Build a legend for graphs appearing later.
RobustIndex = 1;
MixedIndex = 1;
Legend = cell(NumOfAlgorithms,1);
for i=1:length(Legend)
    switch AllScenariosData(1).Data(i).AlgorithmType
        case -1
            Legend{i} = 'Benchmark';
        case 0
            Legend{i} = 'Nominal';
        case 1
            Legend{i} = ['Rob. L_\infty #',num2str(RobustIndex)];
            RobustIndex = RobustIndex + 1;
        case 2
            Legend{i} = ['Rob. Mixed #',num2str(MixedIndex)];
            MixedIndex = MixedIndex + 1;
    end
end

% Build a structure to store statistics
Statistics.PerDay(NumOfScenarios).TrueCost = [];
Statistics.PerDay(NumOfScenarios).TrueOverEstimateCost = [];
if(1==WasBenchmarkRun)
    Statistics.PerDay(NumOfScenarios).CostOverOptimal = [];
end
if(1==WasNominalRun)
    Statistics.PerDay(NumOfScenarios).DifferenceFromNominal = [];
end

Statistics.PerYear(NumOfScenarios).TrueCost = [];
Statistics.PerYear(NumOfScenarios).TrueOverEstimateCost = [];
if(1==WasBenchmarkRun)
    Statistics.PerYear(NumOfScenarios).CostOverOptimal = [];
end
if(1==WasNominalRun)
    Statistics.PerYear(NumOfScenarios).DifferenceFromNominal = [];
end



for i=1:NumOfScenarios
    CurrentScenarioData = AllScenariosData(i).Data;
    TitleSuffix = [' - Scenario #',num2str(i)];
    
    
    
    %Plot "true" cost of schedule - cost of schedule for ground-truth
    %demand
    figure(100*(i-1)+1);
    for j=1:length(CurrentScenarioData)
        plot(CurrentScenarioData(j).TrueCost);
        hold all;
    end
    grid on;
    xlabel('Day in year (from January 15^{th}');
    ylabel('Cost [$]');
    title(['Cost for Ground Truth Demand',TitleSuffix]);
    legend(Legend);
    
    
    %Plot Ratio of Estimated cost over True Cost (Budget exceeding)
    figure(100*(i-1)+2);
    for j=2:length(CurrentScenarioData)
        plot((CurrentScenarioData(j).TrueCost./CurrentScenarioData(j).EstimatedCost-1)*100);
        hold all;
    end
    grid on;
    xlabel('Day in year (from January 15^{th}');
    ylabel('Cost over Expected [%]');
    title(['Cost over Estimate ("Over Budget")',TitleSuffix]);
    legend(Legend(2:end));
    
    for j=1:length(CurrentScenarioData)    
        TrueCost = CurrentScenarioData(j).TrueCost;
        
        Statistics.PerDay(i).TrueCost = TrueCost;
        Statistics.PerYear(i).TrueCost(j).Average = mean(TrueCost);
        Statistics.PerYear(i).TrueCost(j).Std = std(TrueCost);
        Statistics.PerYear(i).TrueCost(j).Max = max(TrueCost);
        Statistics.PerYear(i).TrueCost(j).Min = min(TrueCost);
        Statistics.PerYear(i).TrueCost(j).Median = median(TrueCost);
               
        Ratio = (CurrentScenarioData(j).TrueCost./CurrentScenarioData(j).EstimatedCost-1)*100;
        Statistics.PerDay(i).TrueOverEstimateCost = Ratio;
        Statistics.PerYear(i).TrueOverEstimateCost(j).TimePositive = sum(Ratio > 0)/length(Ratio) * 100;
        Statistics.PerYear(i).TrueOverEstimateCost(j).Max = max(Ratio);
        Statistics.PerYear(i).TrueOverEstimateCost(j).ConditionalValueAtRisk = mean(Ratio(Ratio > 0));
        Statistics.PerYear(i).TrueOverEstimateCost(j).StdofValueAtRisk = std(Ratio(Ratio > 0));
    end
    
    if(1==WasBenchmarkRun)
        BenchmarkData = CurrentScenarioData(1);
        
        %Compute True cost over Nominal True Cost.
        figure(100*(i-1)+3);
        for j=2:length(CurrentScenarioData)
            plot((CurrentScenarioData(j).TrueCost./BenchmarkData.TrueCost-1)*100);
            hold all;
        end
        grid on;
        xlabel('Day in year (from January 15^{th}');
        ylabel('Cost over Benchmark [%]');
        title(['Cost over Benchmark for Ground Truth Demand',TitleSuffix]);
        legend(Legend(2:end));
        
        for j=1:length(CurrentScenarioData)     
            Ratio = (CurrentScenarioData(j).TrueCost./BenchmarkData.TrueCost-1)*100;    
            Statistics.PerDay(i).CostOverOptimal = Ratio;
            if(j==1) %Benchmark algorithm
                Statistics.PerYear(i).CostOverOptimal(j).TimePositive = 0;
                Statistics.PerYear(i).CostOverOptimal(j).Max = 0;
                Statistics.PerYear(i).CostOverOptimal(j).ConditionalValueAtRisk = 0;
                Statistics.PerYear(i).CostOverOptimal(j).StdofValueAtRisk = 0;
            else %Other algorithms
                Statistics.PerYear(i).CostOverOptimal(j).TimePositive = sum(Ratio > 0)/length(Ratio) * 100;
                Statistics.PerYear(i).CostOverOptimal(j).Max = max(Ratio);
                Statistics.PerYear(i).CostOverOptimal(j).ConditionalValueAtRisk = mean(Ratio(Ratio > 0));
                Statistics.PerYear(i).CostOverOptimal(j).StdofValueAtRisk = std(Ratio(Ratio > 0));
            end
        end
    end
    
    
    if(1 == WasNominalRun)
        NominalData = CurrentScenarioData(2);

        %Compute how different the solution is different from the nominal solution.
        figure(100*(i-1)+4);
        for j=3:length(CurrentScenarioData)
            plot(sum((CurrentScenarioData(j).Power_Generation ~= NominalData.Power_Generation |...
                CurrentScenarioData(j).Heat_Generation  ~= NominalData.Heat_Generation  |...
                CurrentScenarioData(j).Fuel_Consumption ~= NominalData.Fuel_Consumption),2));
            hold all;
        end
        grid on;
        xlabel('Day in year (from January 15^{th}');
        ylabel('#hours with different schedule from nominal');
        title(['Difference from Nominal Solution',TitleSuffix]);
        legend(Legend(3:end));
        
        for j=1:length(CurrentScenarioData)

            Diff = sum((CurrentScenarioData(j).Power_Generation ~= NominalData.Power_Generation |...
                CurrentScenarioData(j).Heat_Generation  ~= NominalData.Heat_Generation  |...
                CurrentScenarioData(j).Fuel_Consumption ~= NominalData.Fuel_Consumption),2);
            
            Statistics.PerDay(i).DifferenceFromNominal = Diff;
            Statistics.PerYear(i).DifferenceFromNominal(j).Max = max(Diff);
            Statistics.PerYear(i).DifferenceFromNominal(j).Min = min(Diff);
            Statistics.PerYear(i).DifferenceFromNominal(j).Average = mean(Diff);
            Statistics.PerYear(i).DifferenceFromNominal(j).Median = median(Diff);
        end
    end
end

end