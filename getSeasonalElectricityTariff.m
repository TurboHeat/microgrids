function [dailyTariffs] = getSeasonalElectricityTariff(buildingType, currentDate, dt)
%getSeasonalElectricityTariff Returns a vector with the electricity
%tariff for each time step during the considered day, for the requested building. 
% Based on the createElectricityTariffProfile.m file.
%% Handling inputs
arguments
  buildingType (1,1) BuildingType {mustBeInteger, mustBeGreaterThanOrEqual(buildingType,1), mustBeLessThanOrEqual(buildingType,5)}
  currentDate (1,1) datetime
  dt (1,1) double {mustBeInteger, mustBePositive} = 15 % length of time step in [s]
end
if buildingType == 5
  buildingType = 1; % If the building is hospital use "commercial tall", which is
end                 % the same tariff as large hotel (buildingType=1)

%% Constants
SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
HOURS_PER_DAY = 24;
RATES = LOAD_TARIFFS();

%% Actual computation...
stepsPerHour = SECONDS_PER_HOUR / dt; %number of lines per hour
tariffQueryTimes = hours(linspace(0, HOURS_PER_DAY, stepsPerHour*HOURS_PER_DAY+1).'); tariffQueryTimes(end) = [];

switch buildingType
  case BuildingType.LargeHotel    
    if day(currentDate) ~= Weekday.Sunday && any(month(currentDate) == 6:9)
      rate = RATES("285S");
    else
      rate = RATES("285");
    end
    
  case {BuildingType.FullServiceRestaurant, BuildingType.SmallHotel} % aka "FSR"
    if all(day(currentDate) ~= [Weekday.Sunday, Weekday.Saturday]) && any(month(currentDate) == 6:9)
      rate = RATES("282S");
    else
      rate = RATES("282");
    end

  case BuildingType.ResidentialHIGH % Residential high
    if any(day(currentDate) == [Weekday.Sunday, Weekday.Saturday])
      rate = RATES("184WE");
    else
      if any(month(currentDate) == 6:9)
        rate = RATES("184S");
      else
        rate = RATES("184W");
      end
    end
end
dailyTariffs = rate.tariffAtTime(tariffQueryTimes);
end

function tariffs = LOAD_TARIFFS()
% All data is taken from https://www.lipower.org/about-us/tariff/
% Last update: LIPA-Tariff-September-2020.pdf
%% Rate 182 (unused)

%% Rate 184 (June-September Weekdays)
% Voluntary Large Residential Service with Multiple Rate Periods
% Here we are ignoring the "first 125 kWh" rates!
OOI = [0.0256, 0.2853, NaN]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(10,1); ...  0am-10am
                2*ones(10,1); ... 10am- 8pm
                  ones(4,1)];  %   8pm- 0am
R184S = ElectricityTariff( OOI(tariffIdAtHour), "184S", (0:23).' );

%% Rate 184 (October-May Weekdays)
% Voluntary Large Residential Service with Multiple Rate Periods
% Here we are ignoring the "first 125 kWh" rates!
OOI = [0.0256, 0.0801, NaN]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(10,1); ...  0am-10am
                2*ones(10,1); ... 10am- 8pm
                  ones(4,1)];  %   8pm- 0am
R184W = ElectricityTariff( OOI(tariffIdAtHour), "184W", (0:23).' );

%% Rate 184 (Weekends)
% Voluntary Large Residential Service with Multiple Rate Periods
OOI = [0.0256, NaN, NaN]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = ones(24,1);
R184WE = ElectricityTariff( OOI(tariffIdAtHour), "184WE", (0:23).' );

%% Rate 282 (June-September Weekdays)
% Voluntary Large Demand Metered Service With Multiple Rate Periods
OOI = [0.0035, 0.0251, 0.021]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(7,1); ...  0am- 7am
                3*ones(5,1); ...  7am-12pm
                2*ones(8,1); ... 12pm- 8pm
                3*ones(3,1); ...  8pm-11pm
                  ones(1,1)]; %  11pm- 0am
R282S = ElectricityTariff( OOI(tariffIdAtHour), "282S", (0:23).' );

%% Rate 282 (Non-Summer + Weekends in Summer)
% Voluntary Large Demand Metered Service With Multiple Rate Periods
OOI = [0.0035, NaN, 0.021]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(7,1);  ...  0am- 7am
                3*ones(16,1); ...  7am-11pm
                  ones(1,1)];  %  11pm- 0am
R282 = ElectricityTariff( OOI(tariffIdAtHour), "282", (0:23).' );

%% Rate 285 (June-September, "Secondary" consumer)
% Large General and Industrial Service With Multiple Rate Periods
OOI = [0.0058, 0.0376, 0.024]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(7,1);  ...  0am- 7am
                3*ones(3,1);  ...  7am-10am
                2*ones(12,1); ... 10am-10pm
                3*ones(2,1)]; %   10pm- 0am
R285S = ElectricityTariff( OOI(tariffIdAtHour), "285 (Summer)", (0:23).' );

%% Rate 285 (October-May + Sunday in Summer)
% Large General and Industrial Service With Multiple Rate Periods
OOI = [0.0058, NaN, 0.024]; % Off-Peak/On-Peak/Intermediate rates
tariffIdAtHour = [ones(7,1); ...  0am-7am
                3*ones(17,1)]; %  7am-0am
R285 = ElectricityTariff( OOI(tariffIdAtHour), "285", (0:23).' );

%% Rate 580 (unused)

%% Pack all tariffs
tariffs = containers.Map(["184S", "184W", "184WE", "282S", "282", "285S", "285"],...
                         {R184S,  R184W,  R184WE,  R282S,  R282,  R285S,  R285});
end
