function el_power_plot(power_map)

valve_pos = 0:0.2:0.8;
% speed = []
x = [];
speed_level = [];
z = power_map(1:end-1);
for i = 1:length(power_map) - 1
  valve_pos_level = floor((i -1)/9) + 1;
  x = [x, valve_pos(valve_pos_level)]
  speed_level = [speed_level, mod(i-1, 9) + 1]
end

plot3(speed_level, x, z, '.'), xlabel('speed level'), ylabel('valve pos'), zlabel('el. power')