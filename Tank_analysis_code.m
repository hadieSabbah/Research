%% New code %%

%% Tank discharge Code %%

%Variables
P_0_tank = 1.379e+6; %pascals 
P_f= 448159; %pascals
T_0 = 300; %kelvin
P_list = [P_0_tank/P_f:-0.1:1];
gamma = 1.4;
R = 287; %j/kg*k
r = 1/2;%m
h = 4;
A_star_old = 0.0088; 
V_tank_new = pi*r^2*h;
V_tank_old = 10.3;

%Prelocating Variables 
M = [];
P_tank_cur = [];
t_old = [];
C_d = 0.995;
for i = 1: numel(P_list)
 [M(i),T,P,rho,area] = flowisentropic(gamma,P_list(i)); 
 a_0 = ((gamma*R*T_0)/M(i))^(1/2);
 tau_old = (V_tank_old/(A_star_old*a_0))*((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));
 t_old(i) = log((1/P_list(i)))*(-tau_old);
 P_tank_cur(i) = (1/P_list(i))*P_0_tank;
 rho_tank_cur(i) = P_tank_cur(i)/(T_0*R);
 m_flow(i) = C_d*A_star_old*rho_tank_cur(i)*sqrt(((gamma*R*T_0)/M(i)))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
end

t_final = t_old(end);

figure(1)
plot(t_old,(P_tank_cur));
grid on 
xlabel("Time")
ylabel("$P_{tank}$","Interpreter","latex");
title("Tank Pressure vs Time");

%Mass flow rate vs Time