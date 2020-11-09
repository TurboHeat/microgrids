function [Output] = RunNominalAlgorithm(timeStepSize, endTime, kwargs)
% Analysis of the solution found by the shortestpath solver for no
% uncertainty, appearing on (Rist et al., 2017).
%
arguments
  timeStepSize (1,1) double {mustBePositive} = 15 % length of a timestep in [s]
  endTime (1,1) double {mustBePositive} = 24; % final time in [h]
                            
  % Parameter for first uncertainty set. 
  % In equation (6) in the paper, we take:
  %   Delta_P(t) = alpha_Linfty*std_P(t) 
  %   Delta_H(t) = alpha_Linfty*std_H(t).
  
  kwargs.PriceIndex (1,1) double {mustBePositive} = 2;
  kwargs.BuildingType (1,1) double {mustBePositive} = 2;
  
    
  kwargs.dataPath (1,1) string = "../Data"
  kwargs.transitionPenalty (1,1) double = 0.01;
end
PriceIndex = kwargs.PriceIndex;
BuildingType = kwargs.BuildingType;
dataPath = kwargs.dataPath;
transitionPenalty = kwargs.transitionPenalty;
Output = RunRobustLinftyAlgorithm(timeStepSize,endTime,'alpha',0,...
                                  'dataPath',dataPath,...
                                  'transitionPenalty',transitionPenalty,...
                                  'PriceIndex',PriceIndex,...
                                  'BuildingType',BuildingType);

Output.AlgorithmType = 0; 
Output.AlgorithmParameters{1} = [];
end