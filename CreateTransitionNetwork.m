function CreateTransitionNetwork(kwargs)
%--------------------------------------------------------------%
% TODO
%--------------------------------------------------------------%

%% Handle inputs:
arguments % accepted key-value pairs:
  kwargs.smin (1,1) uint8 = 1;  % minimum 'on' speed level
  kwargs.smax (1,1) uint8 = 25; % maximum speed level
  kwargs.vmin (1,1) uint8 = 0;  % minimum bypass valve level (corresponds to 0% )
  kwargs.vmax (1,1) uint8 = 9;  % maximum bypass valve level (corresponds to 90%)
  kwargs.dt   (1,1) double = 15; % duration of a step [s]
  kwargs.endTime (1,1) uint16 = 24*60; % final time [min]. Default: 24 [h/day] * 60 [min/h]  
  kwargs.savePath (1,1) string = fullfile(fileparts(mfilename), "Data");
end

tic
% Unpack inputs:
smin = kwargs.smin;
smax = kwargs.smax;
vmin = kwargs.vmin;
vmax = kwargs.vmax;
dt = kwargs.dt;
endTime = kwargs.endTime;
savePath = kwargs.savePath;

%% Constants:
SECONDS_PER_MINUTE = 60;

%% Preliminary computations:
tStartup = 120/dt + 8 * 2; % [timesteps]
tShutdown = 180/dt; % [timesteps]
stepsPerMinute = SECONDS_PER_MINUTE / dt;
nt = single(endTime * stepsPerMinute); % number of time steps
ns = (smax - smin + 1);
nv = (vmax - vmin + 1);
nWorkingStates = ns * nv; % number of possible nodes per time step, excluding the 'off' state
nStates = single(nWorkingStates + 1);
sm = StateMapper32( ceil(log2(nStates)), ceil(log2(nt)) );

%% Generate all states at all times:
[SS,VV] = ndgrid(uint16(smin:smax), uint16(vmin:vmax));
svToStateNumber = [0,0; [SS(:), VV(:)]]; % this will be used in the uint32 mapping

[SS,TT] = ndgrid(1:nStates, 1:nt);
tSr = [zeros(1,3,'uint16');                 % Special row indicating "start"
       intmax('uint16')*ones(1,3,'uint16'); % Special row indicating "finish"
       [TT(:), SS(:), nt-TT(:)]]; % formerly, "aux_table"
toc
error('WORK IN PROGRESS!');

%% Create transitions matrix
% Preallocation:
nNodes = size(tSr, 1); %#ok<UNRCH>
ubTransitionsPerTimestep = nWorkingStates * 11; % Need to replace 11 with a f(ns,nv,...)
connectivity = spalloc(nNodes, nNodes, ubTransitionsPerTimestep);

%{
Allowed transitions:
•	(s,v) -> (s,v) - keep same state, advances time by 1.
•	(s,v) -> (s,v+1) - increase bypass valve position by 1, advances time by 1. Available only if v<=v_max - 1.
•	(s,v) -> (s,v-1) - decrease bypass valve position by 1, advances time by 1. Available only if v>=v_min + 1.
•	(s,v) -> (s+1,v) - increase engine velocity by 1, advance time by 2. Available only if s<=s_max - 1.
•	(s,v) -> (s-1,v) - decrease engine velocity by 1, advance time by 1. Available only if s >= s_min + 1.
•	(s,v) -> (s+1,v+1) - increase engine velocity by 1 and bypass position by 1, advance time by 2. Available only if s<=s_max - 1 and v<=v_max - 1.
•	(s,v) -> (s+1,v+2) - increase engine velocity by 1 and bypass position by 2, advance time by 2. Available only if s<=s_max - 1 and v<=v_max - 2.
•	(s,v) -> (s+1,v-1) - increase engine velocity by 1 and decrease bypass position by 1, advance time by 2. Available only if s<=s_max - 1 and v>= v_min - 1.
•	(s,v) -> (s+1,v-2) - increase engine velocity by 1 and decrease bypass position by 2, advance time by 2. Available only if s<=s_max - 1 and v>= v_min - 2.
•	(s,v) -> (s-1,v+1) - decrease engine velocity by 1 and increase bypass position by 1, advance time by 1. Available only if s>=s_min- 1 and v<=v_max - 1.
•	(s,v) -> (s-1,v-1) - decrease engine velocity by 1 and bypass position by 1, advance time by 1. Available only if s>=s_min- 1 and v>=v_min+1.
•	off -> (s_max,v) - start up, advance time by T_startup. The bypass valve position can be anything between v_min and v_max.
•	(s_min,v_min) -> off - shutdown, advance time by T_shutdown.
%}

% Transition table info:
% State @ t (state_t); State @ t+1 (state_t1);  Remaining time (RT)

% Transitions from start: any of the 46 possible states
node_table = table(); % TEMPORARY, to avoid a "compilation-time" error.
Transition.state_t = repmat({'Start'}, (nWorkingStates + 1), 1);
Transition.state_t1 = node_table.Name(3:(nWorkingStates + 3)); %%% why does 3 come here?
Transition.cost = zeros(nWorkingStates+1, 1);

Transition_t = struct2table(Transition);

t_varnames = {'state_t', 'state_t1', 'cost'};
%Transitions from any state

tic
parfor jj = 3:size(node_table, 1) % start and end nodes are treated separately
  current_node = node_table(jj, :);
  isoff = ~(isempty(cell2mat(strfind(current_node.Name, 'x')))); %1 when state is off
  if current_node.RT >= 1 %Base condition to change state, otherwise go to end state
    %keep same state, only increase time step
    %this if is to distinguish off states from regular states
    if isoff %1 when state is off
      next_name = {[num2str(current_node.t+1), 'x']};
    else
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s), 'V', num2str(current_node.v)]};
    end
    
    new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
    
    %decrease speed (RT>=1, s>smin)
    if (current_node.s > smin)
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s-1), 'V', num2str(current_node.v)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
    %increase bypass (RT>=1, v<vmax, state=on)
    if (current_node.v < vmax && isoff == 0)
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s), 'V', num2str(current_node.v+1)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
    %decrease bypass (RT>=1, v>vmin)
    if (current_node.v > vmin)
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s), 'V', num2str(current_node.v-1)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
    %decrease speed and increase bypass (RT>=1, s>smin, v<vmax, state=on)
    if (current_node.s > smin && isoff == 0 && current_node.v < vmax)
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s-1), 'V', num2str(current_node.v+1)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
    %decrease speed and decrease bypass (RT>=1, s>smin, v>vmin)
    if (current_node.s > smin && current_node.v > vmin)
      next_name = {['t', num2str(current_node.t+1), 'S', num2str(current_node.s-1), 'V', num2str(current_node.v-1)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
  else % if RT=0, then connect to end state
    next_name = {'End'};
    end_transition = table(current_node.Name, next_name, 0, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; end_transition];
  end
  %increase speed, requires 2 time steps
  if (current_node.RT >= 2 && current_node.s < smax && isoff == 0)
    next_name = {['t', num2str(current_node.t+2), 'S', num2str(current_node.s+1), 'V', num2str(current_node.v)]};
    new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
  end
  %increase both speed (by 1) and bypass by 1 or 2 (RT>=2, s<smax, v<vmax)
  if (current_node.RT >= 2 && current_node.s < smax && isoff == 0 && current_node.v < vmax) %condition to increase bypass 1 unit
    next_name = {['t', num2str(current_node.t+2), 'S', num2str(current_node.s+1), 'V', num2str(current_node.v+1)]};
    new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
    if current_node.v + 1 < vmax %in addition, check if we can increase bypass by 2
      next_name = {['t', num2str(current_node.t+2), 'S', num2str(current_node.s+1), 'V', num2str(current_node.v+2)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
  end
  %increase speed (by 1) and decrease bypass by 1 or 2 (RT>=2, s<smax, v>vmin)
  if (current_node.RT >= 2 && current_node.s < smax && isoff == 0 && current_node.v > vmin) %condition to increase bypass 1 unit
    next_name = {['t', num2str(current_node.t+2), 'S', num2str(current_node.s+1), 'V', num2str(current_node.v-1)]};
    new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
    if current_node.v - 1 > vmin %in addition, check if we can decrease bypass by 2
      next_name = {['t', num2str(current_node.t+2), 'S', num2str(current_node.s+1), 'V', num2str(current_node.v-2)]};
      new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
      Transition_t = [Transition_t; new_transition];
    end
  end
  %startup procedure
  if (isoff == 1 && current_node.RT >= T_startup / dt)
    next_name = {[repmat('t', (vmax + 1), 1), repmat(num2str(current_node.t+T_startup/dt), ...
      (vmax + 1), 1), repmat('S', (vmax + 1), 1), repmat(num2str(smax), (vmax + 1), 1), repmat('V', (vmax + 1), 1), num2str([vmin:vmax]')]};
    new_transition = table(repmat(current_node.Name, (vmax + 1), 1), cellstr(next_name{1, 1}), Inf((vmax + 1), 1), 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
  end
  
  %shutdown procedure
  if (current_node.s == 1 && current_node.v == 0 && current_node.RT >= T_shutdown / dt)
    next_name = {[num2str(current_node.t+T_shutdown/dt), 'x']};
    new_transition = table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
    Transition_t = [Transition_t; new_transition];
  end
end
toc
%
g = digraph(Transition_t.state_t, Transition_t.state_t1, Transition_t.cost);
%
% figure(1);
% h=plot(g);
% layout(h,'layered','Direction','right','Source',[{'Start'}], 'Sink', [{'End'}]);
% highlight(h, [{'Start'}, {'End'}], 'NodeColor', 'r', 'MarkerSize', 6)
% labelnode(h, [{'Start'}, {'End'}], {'Start', 'End'})

%%
savename = ['graph', num2str(endTime/60), 'h.mat'];
save([savePath, savename], 'g', 'node_table', 'Transition_t', 'total_nodes');
save([savePath, 'workspace', '_', savename]);