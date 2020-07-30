function g = CreateTransitionNetwork(kwargs) 
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
tStartup = 120 / dt + 8 * 2; % [timesteps]
tShutdown = 180 / dt; % [timesteps]
stepsPerMinute = SECONDS_PER_MINUTE / dt; 
nt = single(endTime*stepsPerMinute); % number of time steps
ns = (smax - smin + 1); 
nv = (vmax - vmin + 1); 
nWorkingStates = ns * nv; % number of possible nodes per time step, excluding the 'off' state
nStates = single(nWorkingStates+1); 

%% Generate all states
[SS, VV] = ndgrid(uint16(smin:smax), uint16(vmin:vmax)); 
svToStateNumber = [0, 0; [SS(:), VV(:)]]; % this will be used in the uint32 mapping

%% Create adjacency (transitions) matrix
B = BuildStateAdjacencyMatrix(svToStateNumber, tStartup, tShutdown); 
% figure(); spy(reshape(B, nStates, []));

% Compute the indices:
[x, y] = find(reshape(B, nStates, [])); 
x = reshape(x + (0:nStates:nStates * (nt - 1)), [], 1);
y = reshape(y + (nStates:nStates:nStates * nt), [], 1);
A = sparse(x, y, 1); % "1" can be replaced with "true" if we want to save space
A = A(:, 1:size(A,1));
% figure(); spy(B);

%% Compute costs
% TODO!!!

%% Encode transitions as uint32:
% TODO: Figure out if we need this step
%{
sm = StateMapper32(ceil(log2(nStates)), ceil(log2(nt)));
[fromIdx, toIdx] = find(A);
fromTime = ceil(fromIdx / nStates);
% toTime = ceil(toIdx / nStates);
fromState = mod(fromIdx-1, nStates)+1;
toState = mod(toIdx-1, nStates)+1;
uTransitions = sm.toUint(fromState, toState, fromTime);
%}
%{
[SS, TT] = ndgrid(1:nStates, 1:nt); 
tSr = [zeros(1, 3, 'uint16'); ... % Special row indicating "start"
[TT(:), SS(:), nt - TT(:)]; ... 
intmax('uint16') * ones(1, 3, 'uint16')];% Special row indicating "finish"
%{
tSr = [nan(2,3); [TT(:), SS(:), nt-TT(:)]]; % formerly, "aux_table"
%}
%}

%% Create di(rected )graph
g = digraph(A); 
%{
hF = figure(); hAx = axes(hF); hG = plot(hAx, g);
layout(h,'layered','Direction','right','Source',[{'Start'}], 'Sink', [{'End'}]);
highlight(h, [{'Start'}, {'End'}], 'NodeColor', 'r', 'MarkerSize', 6)
labelnode(h, [{'Start'}, {'End'}], {'Start', 'End'})
%}

%% Save results
savename = ['graph_', num2str(endTime/60), 'h.mat']; 
save([savePath, savename], 'g', 'A'); 