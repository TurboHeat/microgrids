mass_flow_rate = [m_dot_f_map(1, :), m_dot_f_map(2, :), m_dot_f_map(3, :), m_dot_f_map(4, :), m_dot_f_map(5, :)]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(mass_flow_rate, power_map(1:end-1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(mass_flow_rate, heat_map(1:end-1))