%--------------------------------------------------------------%
% File: CreateElectricityTariffProfile.m (function)
% Author: Miguel Dias
% Date 14/08/16
% v1.0
% Description: Returns a vector with the electricity tariff for each time
% step during the considered day, for the requested building. Building and
% day options:
% Buildingtypes:
%           1. Large Hotel (commercial tall)
%           2. Full Service Restaurant (commercial medium)
%           3. Small Hotel (commercial medium)
%           4. Residential High (residential rate)
%           5. Hospital (commercial tall)
% Days:     1. - 10.th January - winter 
%           2. - 10.th April - winter
%           3. - 10.th July - summer
%--------------------------------------------------------------%
function [ tariff_vec ] = CreateElectricityTariffProfile( building_type,day,dt)
%CreateElectricityTariffProfile Returns a vector with the electricity
%tariff for each time step during the considered day, for the
% requested building. ONLY FPE, i.e. A parcel off electricity charge. B and
% C costs are added a posteriori
n_lines=3600/dt; %number of lines per hour
tariff_vec=ones(24*n_lines,1);

if building_type==5
    building_type=1; % If the building is hospital, commercial tall,
end                  % apply same tariff as large hotel so bt=1

switch building_type
    case 1 %Large Hotel
       % Commercial tall rate, 285
       %Peak 10am-10pm during the summer, except sunday
       %Off_peak 12pm - 7am
       %Intermediate 7am-10am && 10pm-12pm
       if day==3 %summer
           peak_h=n_lines.*[10 22]; %start and end indices of peak hours
           off_peak_h=[1 n_lines*7]; % one is the start index, midnight
           off_peak_charge     = 0.0291; %$/kWh
           peak_charge         = 0.0543; %$/kWh
           intermediate_charge = 0.0435; %$/kWh
           tariff_vec=intermediate_charge.*tariff_vec;
           tariff_vec([peak_h(1):peak_h(2)])=tariff_vec([peak_h(1):peak_h(2)]).*peak_charge/intermediate_charge; %assign peak charge
           tariff_vec([off_peak_h(1):off_peak_h(2)])=tariff_vec([off_peak_h(1):off_peak_h(2)]).*off_peak_charge/intermediate_charge; %assign off peak charge 
       else %winter
          %No peak hour - only off peak and intermediate
           off_peak_h=[1 n_lines*7]; % one is the start index, midnight
           off_peak_charge     = 0.0291; %$/kWh
           intermediate_charge = 0.0435; %$/kWh
           tariff_vec=intermediate_charge.*tariff_vec;
           tariff_vec([off_peak_h(1):off_peak_h(2)])=tariff_vec([off_peak_h(1):off_peak_h(2)]).*off_peak_charge/intermediate_charge; %assign off peak charge
       end
       
    case 2  %FSR
       % Commercial medium rate, 282
       %Peak 12am-8pm during the summer weekdays
       %Off_peak 11pm - 7am
       %Intermediate 7am-12am && 8pm-11pm
       if day==3 %summer
           peak_h=n_lines.*[12 20]; %start and end indices of peak hours
           off_peak_h=[1 n_lines*7 n_lines*23 n_lines*24]; % one is the start index, midnight
           off_peak_charge     = 0.0273; %$/kWh
           peak_charge         = 0.0444; %$/kWh
           intermediate_charge = 0.0412; %$/kWh
           tariff_vec=intermediate_charge.*tariff_vec;
           tariff_vec([peak_h(1):peak_h(2)])=tariff_vec([peak_h(1):peak_h(2)]).*peak_charge/intermediate_charge; %assign peak charge
           tariff_vec([off_peak_h(1):off_peak_h(2)])=tariff_vec([off_peak_h(1):off_peak_h(2)]).*off_peak_charge/intermediate_charge; %assign 1st off peak period charge 
           tariff_vec([off_peak_h(3):off_peak_h(4)])=tariff_vec([off_peak_h(3):off_peak_h(4)]).*off_peak_charge/intermediate_charge; %assign 2nd off peak period charge
       else %winter
          %No peak hour - only off peak and intermediate
           int_h=n_lines.*[7 23]; % intermediate hours in winter comprise the 'peak' and intermediate summer hours
           off_peak_charge     = 0.0273; %$/kWh
           intermediate_charge = 0.0412; %$/kWh
           tariff_vec=off_peak_charge.*tariff_vec; %assing off peak charge
           tariff_vec([int_h(1):int_h(2)])=tariff_vec([int_h(1):int_h(2)]).*intermediate_charge/off_peak_charge; %assign intermediate charge
       end
    case 3 % Small Hotel - same rate as FSR
               % Commercial medium rate, 282
       %Peak 12am-8pm during the summer weekdays
       %Off_peak 11pm - 7am
       %Intermediate 7am-12am && 8pm-11pm
       if day==3 %summer
           peak_h=n_lines.*[12 20]; %start and end indices of peak hours
           off_peak_h=[1 n_lines*7 n_lines*23 n_lines*24]; % one is the start index, midnight
           off_peak_charge     = 0.0273; %$/kWh
           peak_charge         = 0.0444; %$/kWh
           intermediate_charge = 0.0412; %$/kWh
           tariff_vec=intermediate_charge.*tariff_vec;
           tariff_vec([peak_h(1):peak_h(2)])=tariff_vec([peak_h(1):peak_h(2)]).*peak_charge/intermediate_charge; %assign peak charge
           tariff_vec([off_peak_h(1):off_peak_h(2)])=tariff_vec([off_peak_h(1):off_peak_h(2)]).*off_peak_charge/intermediate_charge; %assign 1st off peak period charge 
           tariff_vec([off_peak_h(3):off_peak_h(4)])=tariff_vec([off_peak_h(3):off_peak_h(4)]).*off_peak_charge/intermediate_charge; %assign 2nd off peak period charge
       else %winter
          %No peak hour - only off peak and intermediate
           int_h=n_lines.*[7 23]; % intermediate hours in winter comprise the 'peak' and intermediate summer hours
           off_peak_charge     = 0.0273; %$/kWh
           intermediate_charge = 0.0412; %$/kWh
           tariff_vec=off_peak_charge.*tariff_vec; %assing off peak charge
           tariff_vec([int_h(1):int_h(2)])=tariff_vec([int_h(1):int_h(2)]).*intermediate_charge/off_peak_charge; %assign intermediate charge
       end
    case 4 % Residential high
        %Determine summer or winter rate - This corresponds to rate 580,
        %the easiest that only distinguishes summer/winter
%         if day == 3 %summer rate, for july
%             charge=0.1071; %$/kWh
%         else % winter rates, for january and april days
%             charge=0.0607; %$/kWh
%         end
%         tariff_vec=charge*tariff_vec; %24h period
        %% Rate 182
        %Peak 10am-8pm
        %Off_peak 8pm - 10am
%         peak_h=n_lines.*[10 20]; %start and end indices of peak hours
%         if day==3 %summer
%             off_peak_charge = 0.0677; %$/kWh
%             peak_charge     = 0.1327; %$/kWh
%         else %winter
%             off_peak_charge = 0.0524; %$/kWh
%             peak_charge     = 0.0524; %$/kWh
%         end
%         tariff_vec=off_peak_charge.*tariff_vec;
%         tariff_vec([peak_h(1):peak_h(2)])=tariff_vec([peak_h(1):peak_h(2)]).*peak_charge/off_peak_charge; %assign peak hour charge
        %% Rate 184
        %Peak 10am-8pm
        %Off_peak 8pm - 10am
        peak_h=n_lines.*[10 20]; %start and end indices of peak hours
        if day==3 %summer
            off_peak_charge = 0.0442; %$/kWh
            peak_charge     = 0.2461; %$/kWh
        else %winter
            off_peak_charge = 0.0442; %$/kWh
            peak_charge     = 0.0866; %$/kWh
        end
        tariff_vec=off_peak_charge.*tariff_vec;
        tariff_vec([peak_h(1):peak_h(2)])=tariff_vec([peak_h(1):peak_h(2)]).*peak_charge/off_peak_charge; %assign peak hour charge
      
end
end

