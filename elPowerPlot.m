function elPowerPlot(power_map)

valve_pos = 0:0.2:0.8;
np = size(power_map,2);
[x, speed_level] = deal(zeros(1,np));
z = power_map(1:end-1);

for i = 1:np - 1
  valve_pos_level = floor((i -1)/9) + 1;
  x(i) = [x, valve_pos(valve_pos_level)];
  speed_level(i) = [speed_level, mod(i-1, 9) + 1];
end

plot3(speed_level, x, z, '.'), xlabel('speed level'), ylabel('valve pos'), zlabel('el. power')