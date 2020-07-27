%--------------------------------------------------------------%
% File: extract_path.m (function)
% Author: Miguel Dias
% Date 15/08/16
% v1.0
% Description: Get MGT commitment from solution of the shortestpath
% algorithm
%[power_MGT, heat_MGT] = extract_path(sol_path)
%--------------------------------------------------------------%
function [power_MGT, heat_MGT, fuel_MGT] = extract_path(path, power_map, heat_map, fuel_map, SV_states)
%EXTRACT_PATH: Get MGT commitment from solution of the shortestpath
% algorithm
%Normalize transition
path = strrep(path, 'Start', 't0S0V0');
path = strrep(path, 'x', 'S0V0');
path = strcat('t', path);
path = strrep(path, 'tt', 't');
path = strrep(path, 'End', [num2str(intmax('uint16')), 'S0V0']);


dT = regexp(path.', '[tSV]', 'split');
dT = vertcat(dT{:});
dT = dT(:, 2:4);
dT = uint16(reshape(sscanf(sprintf('%s*', dT{:}), '%u*'), [], 3)); %{t,s,v}
[~, path] = ismember(uint8(dT(2:size(dT, 1)-1, 2:3)), SV_states, 'rows'); %select all points except 1st and last one (End state)

%path=ceil(nakeinterp1(double(dT(2:end-1,1)),path,(1:double(dT(end-1,1))).')); %make vector the correct lenght and round up to remove 'half states'
%index power map with v,s at each time step
power_MGT = nakeinterp1(double(dT(2:end-1, 1)), power_map(path), (1:double(dT(end-1, 1))).');
heat_MGT = nakeinterp1(double(dT(2:end-1, 1)), heat_map(path), (1:double(dT(end-1, 1))).');
fuel_MGT = nakeinterp1(double(dT(2:end-1, 1)), fuel_map(path), (1:double(dT(end-1, 1))).');
end
