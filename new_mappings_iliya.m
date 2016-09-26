%--------------------------------------------------------------%
% File: assign_costs.m (script)
% Author: Iliya Romm
% Date 15/08/16
% v1.0
% Description: Tranform the Transition table generated by
% CreateTransitionNetwork
% into an equivalent description that takes less space.
% Returns:
% SV_states - maps SV values (columns) to a linear index (in rows).
% time_from - time step at of each 'from' transition
% dT - number of time steps of the transition
% from_state_map - maps number assigned by SV_states to the 'from' node
% from_state_map - maps number assigned by SV_states to the 'to' node
% time_map - maps time step (not used)
%--------------------------------------------------------------%
function [SV_states,time_from,dT,from_state_map,to_state_map] = new_mappings_iliya(Transition_t)
%% Generate mapping:
SV_states = uint8([combvec(1:9,0:4).'; 0 0]);
%% "Normalize" Transition_t:
%Replaces the x of the off state with SOVO and start with t0SOVO 
Transition_t(:,1) = strrep(table2cell(Transition_t(:,1)),'Start','t0S0V0');
Transition_t(:,1:2) = strrep(table2cell(Transition_t(:,1:2)),'x','S0V0');
Transition_t(:,1:2) = strcat('t',table2cell(Transition_t(:,1:2))); Transition_t(:,1:2) = strrep(table2cell(Transition_t(:,1:2)),'tt','t');
Transition_t(:,2) = strrep(table2cell(Transition_t(:,2)),'End',[num2str(intmax('uint16')) 'S0V0']);
%% Parse the "string" states into numbers (happens once!):
%{
% We are using "startTime" as a temporary variable
startTime = regexp(table2cell(Transition_t(:,1:2)),'[tSV]','split'); 
startTime = vertcat(startTime{:}); startTime = startTime(:,2:4); 
startTime = uint16(reshape(sscanf(sprintf('%s*', startTime{:}), '%u*'),[],3));
startTime(:,2) = startTime(:,2);

endTime   = startTime(size(startTime,1)/2+1:end,:);
startTime = startTime(1:size(startTime,1)/2,:);

[~,from_state_map] = ismember(uint8(startTime(:,2:3)),SV_states,'rows'); from_state_map = uint8(from_state_map);
[~,to_state_map]   = ismember(uint8(  endTime(:,2:3)),SV_states,'rows');   to_state_map = uint8(to_state_map);
assert(~any(~from_state_map) && ~any(~to_state_map), 'There shoud be no zero states!')

startTime = startTime(:,1);
endTime   =   endTime(:,1);
% Only now startTime contains what its name implies
%}
% We are using "dT" as a temporary variable
dT = regexp(table2cell(Transition_t(:,1:2)),'[tSV]','split'); 
dT = vertcat(dT{:}); dT = dT(:,2:4); 
dT = uint16(reshape(sscanf(sprintf('%s*', dT{:}), '%u*'),[],3));

[~,from_state_map] = ismember(uint8(dT(1:size(dT,1)/2,2:3)),SV_states,'rows'); 
   from_state_map  = uint8(from_state_map);
[~,to_state_map]   = ismember(uint8(dT(size(dT,1)/2+1:end,2:3)),SV_states,'rows');  
   to_state_map    = uint8(to_state_map);
assert(~any(~from_state_map) && ~any(~to_state_map), 'There shoud be no zero states!')

%time_map = [(uint16(0):dT(size(dT,1)/2)).'; intmax('uint16')];
% Find the appropriate integer class for time:
% if time_map(end-1) < intmax('uint8')
%   time_map(end-1) = uint8(time_map(end-1));
% elseif time_map(end-1) < intmax('uint16')
%   time_map(end-1) = uint16(time_map(end-1));
% elseif time_map(end-1) < intmax('uint32')
%   time_map(end-1) = uint32(time_map(end-1));
% else % n < intmax('uint64')  
%   time_map(end-1) = uint64(time_map(end-1));
% end
time_from = dT(1:size(dT,1)/2,1);
dT = diff(reshape(dT(:,1),[],2),1,2);
% Only now dT contains what its name implies