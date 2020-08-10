%--------------------------------------------------------------%
% Author: Miguel Dias
% Date 14/08/16
% v1.0
% Description:Post processing functions that
% returns the peak demand charges, with and without MGT operation for the
% considered day and building.
% Residential buildings do not have demand charges!!
% Buildingtypes:
%           1. Hospital (commercial tall)
%           2. Full Service Restaurant (commercial medium)
%           3. Small Hotel (commercial medium)
%           4. Residential High (residential rate)
% Days:     1. - 10.th January - winter
%           2. - 10.th April - winter
%           3. - 10.th July - summer
%--------------------------------------------------------------%
function [MGT_PDC, MGT_IDC, ut_PDC, ut_IDC, FC] = generateDemandCharges(d_index, dt, power_demand, new_demand)
tariff_map = uint8([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4; 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]).'; %B, D
building_type = tariff_map(d_index, 1);
day = tariff_map(d_index, 2);
n_lines = 3600 / dt; %number of lines per hour
MGT_PDC = 0; %MGT Peak Demand Charge
MGT_IDC = 0; %MGT Intermediate Demand Charge
ut_PDC = 0; %Utility Peak Demand Charge
ut_IDC = 0; %Utility Intermediate Demand Charge
windowSize = 15 * 60 / dt; %15 min averaging window: 15min*60s/dt
b = (1 / windowSize) * ones(1, windowSize);
a = 1;
power_demand = filter(b, a, power_demand) / 1e3; %filter demand with 15min window averaging, conversion to kW
new_demand = filter(b, a, new_demand) / 1e3;


FC = 0; %Fixed costs=Service charge(per day)+Meter charge(per day)
peak_vec = zeros(24*n_lines, 1);
inter_vec = ones(24*n_lines, 1);
switch building_type
  case 1 %Hospital
    % Commercial tall rate, 285
    %Peak 10am-10pm during the summer, except sunday
    %Off_peak 12pm - 7am
    %Intermediate 7am-10am && 10pm-12pm
    if day == 3 %summer
      peak_h = n_lines .* [10, 22]; %start and end indices of peak hours
      peak_vec(peak_h(1):peak_h(2)) = 1;
      off_peak_h = [1, n_lines * 7]; % one is the start index, midnight
      inter_vec([peak_h(1):peak_h(2), off_peak_h(1):off_peak_h(2)]) = 0;
      peak_charge = 22.44; %$/kW
      intermediate_charge = 5.34; %$/kW
      ut_PDC = max(power_demand.*peak_vec) * peak_charge;
      MGT_PDC = max(new_demand.*peak_vec) * peak_charge;
      ut_IDC = max(power_demand.*inter_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*inter_vec) * intermediate_charge;
      
    else %winter
      %No peak hour - only off peak and intermediate
      ut_PDC = 0;
      MGT_PDC = 0;
      off_peak_h = [1, n_lines * 7]; % one is the start index, midnight
      inter_vec(off_peak_h(1):off_peak_h(2)) = 0;
      intermediate_charge = 5.34; %$/kW
      ut_IDC = max(power_demand.*inter_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*inter_vec) * intermediate_charge;
    end
    FC = 7.66 + 2.5; %$/day
  case 2 %FSR
    % Commercial medium rate, 282
    %Peak 12am-8pm during the summer weekdays
    %Off_peak 11pm - 7am
    %Intermediate 7am-12am && 8pm-11pm
    peak_charge = 45.48; %$/kW
    intermediate_charge = 3.9; %$/kW
    if day == 3 %summer
      peak_h = n_lines .* [12, 20]; %start and end indices of peak hours
      peak_vec(peak_h(1):peak_h(2)) = 1;
      off_peak_h = [1, n_lines * 7, n_lines * 23, n_lines * 24]; % one is the start index, midnight
      inter_vec([peak_h(1):peak_h(2), off_peak_h(1):off_peak_h(2)]) = 0;
      ut_PDC = max(power_demand.*peak_vec) * peak_charge;
      MGT_PDC = max(new_demand.*peak_vec) * peak_charge;
      ut_IDC = max(power_demand.*inter_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*inter_vec) * intermediate_charge;
    else %winter
      %No peak hour - only off peak and intermediate
      ut_PDC = 0;
      MGT_PDC = 0;
      int_h = n_lines .* [7, 23]; % intermediate hours in winter comprise the 'peak' and intermediate summer hours
      peak_vec(int_h(1):int_h(2)) = 1; %easier to construct from vector of zeros, in fact this is inter_vec
      ut_IDC = max(power_demand.*peak_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*peak_vec) * intermediate_charge;
    end
    FC = 1.43 + 0.25; %$/day
  case 3 % Small Hotel - same rate as FSR
    % Commercial medium rate, 282
    %Peak 12am-8pm during the summer weekdays
    %Off_peak 11pm - 7am
    %Intermediate 7am-12am && 8pm-11pm
    peak_charge = 45.48; %$/kW
    intermediate_charge = 3.9; %$/kW
    if day == 3 %summer
      peak_h = n_lines .* [12, 20]; %start and end indices of peak hours
      peak_vec(peak_h(1):peak_h(2)) = 1;
      off_peak_h = [1, n_lines * 7, n_lines * 23, n_lines * 24]; % one is the start index, midnight
      inter_vec([peak_h(1):peak_h(2), off_peak_h(1):off_peak_h(2)]) = 0;
      ut_PDC = max(power_demand.*peak_vec) * peak_charge;
      MGT_PDC = max(new_demand.*peak_vec) * peak_charge;
      ut_IDC = max(power_demand.*inter_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*inter_vec) * intermediate_charge;
    else %winter
      %No peak hour - only off peak and intermediate
      ut_PDC = 0;
      MGT_PDC = 0;
      int_h = n_lines .* [7, 23]; % intermediate hours in winter comprise the 'peak' and intermediate summer hours
      peak_vec(int_h(1):int_h(2)) = 1; %easier to construct from vector of zeros, in fact this is inter_vec
      ut_IDC = max(power_demand.*peak_vec) * intermediate_charge;
      MGT_IDC = max(new_demand.*peak_vec) * intermediate_charge;
    end
    FC = 1.43 + 0.25; %$/day
  case 4 % Residential high
    MGT_PDC = 0; %MGT Peak Demand Charge
    MGT_IDC = 0; %MGT Intermediate Demand Charge
    ut_PDC = 0; %Utility Peak Demand Charge
    ut_IDC = 0; %Utility Intermediate Demand Charge
    FC = 1.65; %$/day
end