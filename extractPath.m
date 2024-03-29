%--------------------------------------------------------------%
% File: extractPath.m (function)
% Author: Miel Sharf
% Date 07/08/2020
% v2.0
% Description: Get MGT commitment from solution of the shortestpath
% algorithm
%[power_MGT, heat_MGT] = extractPath(sol_path)
%--------------------------------------------------------------%
function [power_MGT, heat_MGT, fuel_MGT] = extractPath(path, power_map, heat_map, fuel_map, SV_states)
%extractPath: Get MGT commitment from solution of the shortestpath
% algorithm
nStates = size(SV_states,1);
path_aux = path - 1; %Shift - move source node to index 0.;
path_states = mod(path_aux(2:end-1)-1, nStates);
path_times = floor((path_aux(2:end-1) - 0.5)/nStates) + 1;

path_states(path_states == 0) = size(SV_states, 1);

%path=ceil(nakeinterp1(double(dT(2:end-1,1)),path,(1:double(dT(end-1,1))).')); %make vector the correct lenght and round up to remove 'half states'
%index power map with v,s at each time step
power_MGT = nakeinterp1(double(path_times)', power_map(path_states), (1:double(path_times(end))).');
heat_MGT  = nakeinterp1(double(path_times)', heat_map(path_states),  (1:double(path_times(end))).');
fuel_MGT  = nakeinterp1(double(path_times)', fuel_map(path_states),  (1:double(path_times(end))).');
end
