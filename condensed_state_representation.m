

%--------------------------------------------------------------%
% File: condensed_state_representation.m (function)
% Author: Iliya/Miguel
% Date 18/08/16
% v1.0
% Description: Altrernative description of the state using a uint32,
% where the first 6 bits are the state from, next 6 are the state to and
% remaining bits(MSB, the leftmost ones) represent the timestep.
%--------------------------------------------------------------%
clearvars;
close all;
clc
addpath('C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data');
load graph24h_newSU.mat 'Transition_t';
load graph_data_all_days.mat;

state_from = Transition_t.state_t;
state_to = Transition_t.state_t1;

%% Mapping from time and state to single number
new_from = bitshift(uint32(time_from), 6+6, 'uint32') + uint32(from_state_map); %Time in leftmost 20 bits,
%need to shift it 6+6 bits to accomodate the time to and time from in that
%space. Since the 'from' state bits are in the right, no need to shift
%them. The bits between time and from state will be zero.
new_to = bitshift(uint32(time_from+n_tsteps), 6+6, 'uint32') + bitshift(uint32(to_state_map), 6); %Time in leftmost 20 bits,
%need to shift it 6+6 bits to accomodate the time to and time from in that
%space. The 'to' state bits are right of time, so they need to be shifted
%by 6 bits. The rightmost bits, corresponding to the from state will be 0.

%% Reverse mapping, from number to sv states+time
rev_map_from = bitand(bitshift(new_from, 0), uint32(63)); %63=2^6-1-> and with 6 ones in the right: 111 111
%rev_map_to  = bitand(new_to,uint32(4032)); %4032=2^12-63-1 -> and with 12 ones in the right: 111 111 000 000

rev_map_to = bitand(bitshift(new_to, -6), uint32(63)); %shift 6 to the right and then and with 63

rev_map_time = bitand(bitshift(new_to, -12), uint32(2^20-1));