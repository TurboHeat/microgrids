%--------------------------------------------------------------%
% File: CreateTransitionNetwork.m (script)
% Author: Miguel Dias
% Date 10/08/16
% v4.0
% Description: Final transition network. Updated initial transition
% according to John's email. Changes in line 32 and 139-143.
% 1) Generate table with all possible nodes up to defined time step T
% 2) Append the start and end states to the table
% dt=15s, T=59 time steps (same as John in the example)
% in order to capture full start-up/shutdown procedure
% Include 2 input transitions: can chage s and v in the same step!
% STRUCTURE DEFINITIONS:
% Node=(Name,t,s,v,RT)
% t-timestep; RT - remaining time
% Transition=(state_t, state_t1, cost)
% Includes loop parallelization in the main for loop.
% Does not need any external files.
% Timing stats:
%Elapsed time is 1351.708937 seconds. -> construction of the node table
%Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
%Elapsed time is 8612.214804 seconds. -> construction of edges
%--------------------------------------------------------------%
clearvars; close all; clc;
savepath='C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data\';

tic
end_time= 24*60; %end time in minutes
dt=15; %in s
T=end_time*60/dt; %number of time steps
%T=59;
smin=1; %minimum 'on' speed level
smax=9; %maximum speed level
vmin=0; %minimum bypass valve level (corresponds 0%)
vmax=4; %maximum bypass valve level (corresponds to 80%)

nnodes=smax*(vmax+1); %number of possible nodes per time step, without off state
total_nodes=nnodes+1;
T_startup=120+8*2*dt; %in s
T_shutdown= 180; %in s
varnames={'Name' 't' 's' 'v'};

%generate initial and start nodes
Nodes= table({'Start'; 'End'}, [NaN;NaN],[NaN;NaN], [NaN;NaN],'VariableNames', varnames);
Nodes.RT=[NaN; NaN];

%% Generate 'base' table with all the possible states, for a timestep
svvals=combvec(smin:smax,vmin:vmax)';
names=repmat({'name'},nnodes,1);
basic_table=table(names, ones(nnodes,1), svvals(:,1), svvals(:,2), 'VariableNames', varnames);

%append off state before defining RT and Flags
off= table({'1x'}, 1,0,0,  'VariableNames', varnames);
basic_table=[off; basic_table]; %basic table, with off state and all possible combinations

aux_table=repmat(basic_table, T, 1); %all possible nodes, for all timesteps
%% Define the time step for all nodes
ts=1;
for ii=1:(nnodes+1):size(aux_table,1)
    aux_table.t(ii:size(aux_table,1))=ts;
    ts=ts+1;
end
%% Name all the nodes eg: t1S1V0 -time 1, speed 1, valve postition 0
for jj=1:size(aux_table,1)
    aux_table.Name(jj)={['t' num2str(aux_table.t(jj)) 'S' num2str(aux_table.s(jj)) 'V' num2str(aux_table.v(jj))]};
    if aux_table.s(jj)==0 && aux_table.v(jj)==0
        aux_table.Name(jj)={[num2str(aux_table.t(jj)) 'x']};
    end
end

%define RT  for all nodes
aux_table.RT=T-aux_table.t;

%append start and end nodes - final node table
node_table=[Nodes; aux_table]; %final node table
toc
%% Create Transition table, from nodes

% Allowable transitions:
% increase s, RT>=2 (RT=Remaining Time = T-current time step)
% decrease s, RT>=1
% increase v, RT>=1
% decrease v, RT>=1
% Startup, s=v=0 , RT>=T_startup
% Shutdown, s=1, v=0, RT>=T_shutdown

% Transition table info:
% State @ t (state_t); State @ t+1 (state_t1);  Remaining time (RT)

% Transitions from start: any of the 46 possible states
Transition.state_t=repmat({'Start'}, (nnodes+1),1);
Transition.state_t1=node_table.Name(3:(nnodes+3));
Transition.cost=zeros(nnodes+1,1);

Transition_t=struct2table(Transition);

t_varnames={'state_t', 'state_t1', 'cost'};
%Transitions from any state

tic
parfor jj=3:size(node_table,1) % start and end nodes are treated separately
    current_node= node_table(jj,:);
    isoff=~(isempty(cell2mat(strfind(current_node.Name, 'x')))); %1 when state is off
    if current_node.RT>=1 %Base condition to change state, otherwise go to end state
        %keep same state, only increase time step
        %this if is to distinguish off states from regular states
        if isoff %1 when state is off
            next_name={[num2str(current_node.t+1) 'x']};
        else
            next_name={['t' num2str(current_node.t+1) 'S' num2str(current_node.s) 'V' num2str(current_node.v)]};
        end

        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
        
        %decrease speed (RT>=1, s>smin)
        if (current_node.s>smin)
            next_name={['t' num2str(current_node.t+1) 'S' num2str(current_node.s-1) 'V' num2str(current_node.v)]};
            new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
            Transition_t=[Transition_t; new_transition];
        end
        %increase bypass (RT>=1, v<vmax, state=on)
        if (current_node.v<vmax && isoff==0)
            next_name={['t' num2str(current_node.t+1) 'S' num2str(current_node.s) 'V' num2str(current_node.v+1)]};
            new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
            Transition_t=[Transition_t; new_transition];
        end
        %decrease bypass (RT>=1, v>vmin)
        if (current_node.v>vmin)
            next_name={['t' num2str(current_node.t+1) 'S' num2str(current_node.s) 'V' num2str(current_node.v-1)]};
            new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
            Transition_t=[Transition_t; new_transition];
        end
    else % if RT=0, then connect to end state
        next_name={'End'};
        end_transition=table(current_node.Name, next_name, 0, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; end_transition];
    end
    %increase speed, requires 2 time steps
    if (current_node.RT>=2 && current_node.s<smax && isoff==0)
        next_name={['t' num2str(current_node.t+2) 'S' num2str(current_node.s+1) 'V' num2str(current_node.v)]};
        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
    end
    %increase both speed (by 1) and bypass by 1 or 2 (RT>=2, s<smax, v<vmax)
        if (current_node.RT>=2 && current_node.s<smax && isoff==0 && current_node.v<vmax) %condition to increase bypass 1 unit
        next_name={['t' num2str(current_node.t+2) 'S' num2str(current_node.s+1) 'V' num2str(current_node.v+1)]};
        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
            if current_node.v+1<vmax %in addition, check if we can increase bypass by 2
                next_name={['t' num2str(current_node.t+2) 'S' num2str(current_node.s+1) 'V' num2str(current_node.v+2)]};
                new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
                Transition_t=[Transition_t; new_transition];
            end
        end
        %increase speed (by 1) and decrease bypass by 1 or 2 (RT>=2, s<smax, v>vmin)
        if (current_node.RT>=2 && current_node.s<smax && isoff==0 && current_node.v>vmin) %condition to increase bypass 1 unit
        next_name={['t' num2str(current_node.t+2) 'S' num2str(current_node.s+1) 'V' num2str(current_node.v-1)]};
        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
            if current_node.v-1>vmin %in addition, check if we can decrease bypass by 2
                next_name={['t' num2str(current_node.t+2) 'S' num2str(current_node.s+1) 'V' num2str(current_node.v-2)]};
                new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
                Transition_t=[Transition_t; new_transition];
            end
        end
        %decrease speed (by 1) and increase bypass by 1 or 2 (RT>=2, s<smax, v>vmin)
        if (current_node.RT>=1 && current_node.s>smin && isoff==0 && current_node.v<vmax) %condition to increase bypass 1 unit
        next_name={['t' num2str(current_node.t+1) 'S' num2str(current_node.s-1) 'V' num2str(current_node.v+1)]};
        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
        end
    %startup procedure
    if (isoff==1 && current_node.RT>=T_startup/dt)
        next_name={[repmat('t',(vmax+1),1) repmat(num2str(current_node.t+T_startup/dt),...
        (vmax+1),1) repmat('S',(vmax+1),1) repmat(num2str(smax),(vmax+1),1) repmat('V',(vmax+1),1) num2str([vmin:vmax]')]};
        new_transition=table(repmat(current_node.Name,(vmax+1),1), cellstr(next_name{1,1}), Inf((vmax+1),1), 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
    end
    
    %shutdown procedure
    if (current_node.s==1 && current_node.v==0 && current_node.RT>=T_shutdown/dt)
        next_name={[num2str(current_node.t+T_shutdown/dt) 'x']};
        new_transition=table(current_node.Name, next_name, Inf, 'VariableNames', t_varnames);
        Transition_t=[Transition_t; new_transition];
    end
end
toc
%
g=digraph(Transition_t.state_t, Transition_t.state_t1, Transition_t.cost);
%
% figure(1);
% h=plot(g);
% layout(h,'layered','Direction','right','Source',[{'Start'}], 'Sink', [{'End'}]);
% highlight(h, [{'Start'}, {'End'}], 'NodeColor', 'r', 'MarkerSize', 6)
% labelnode(h, [{'Start'}, {'End'}], {'Start', 'End'})
%%
savename= ['graph', num2str(end_time/60),'h.mat'];
save([savepath,savename], 'g', 'node_table', 'Transition_t', 'total_nodes');
save([savepath, 'workspace','_', savename]);