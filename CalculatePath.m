function [power, heat, cost] = CalculatePath(path, power_map, heat_map, g, SV_states)

cost = 0;
l = length(path);

for i = 1:(l - 1)
  FromNode = path{i};
  ToNode = path{i+1};
  Index = find(ismember(g.Edges.EndNodes(:, 1), FromNode) & ismember(g.Edges.EndNodes(:, 2), ToNode));
  cost = cost + g.Edges.Weight(Index);
end

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
power = nakeinterp1(double(dT(2:end-1, 1)), power_map(path), (1:double(dT(end-1, 1))).');
heat = nakeinterp1(double(dT(2:end-1, 1)), heat_map(path), (1:double(dT(end-1, 1))).');
% fuel_MGT = nakeinterp1(double(dT(2:end-1,1)), fuel_map(path),(1:double(dT(end-1,1))).');
end