%% Calculating mass flow 
P0 = [450:20:690]; %kpa
T0 = 319; %kelvin %Assuming constant stagnation temperature 
mach = 2.5;
R = 287;
gamma = 1.4;
A = 0.02322576;

T_ratio = [];
P_ratio = [];
[mach,T_ratio,P_ratio,rho,area] = flowisentropic(gamma,mach);

T = [];
P = [];
U_inf = [];
m_dot = [];
for i = 1:numel(P0)
    T = T_ratio*T0;
    P(i) = P_ratio*P0(i);
    a_inf = sqrt(gamma*R*T);
    U_inf = mach*a_inf;
    rho(i) = (P(i)*1000)/(R*T);
    m_dot(i) = rho(i)*A*U_inf;
end

%Assuming calorically perfect gas---> gamma = 1.4
plot(P0,m_dot)
xlabel("$P_{0}\:[kPa]$","Interpreter","latex")
ylabel("$\dot{m}\:[Kg/s]$","Interpreter","latex")
title("Mass Flow Vs Stagnation Pressure","Interpreter","latex")
grid on

%Creating table...
m_dot_trans = m_dot.';
P0_trans = P0.';
mass_press_table = table(m_dot_trans,P0_trans);

%% Blowout Math tank Math 

%Initilazing initial equations

%%% ISOTHERMAL CASE %%%%
M = 2.5;
P_tank_max = 1.419e+6; %pascals
P_atm = 101325;
A_star = linspace(0.01,0.1,10);
gamma = 1.4;
T_0 = 300; %kelivn
R = 287;

C_d = 0.995;
V_tank = 10.3;%m^3
a_0 = ((gamma*R*T_0)/M)^(1/2);
tau = [];

for i = 1:numel(A_star)
    tau(i) = (V_tank*((gamma+1)/2)^((gamma+1)/(2*(gamma-1))))/(C_d*A_star(i)*a_0);
end

P_P0 = [0.1:0.1:1];
t_isothermal = [];

for i = 1:numel(A_star)
    t_isothermal(i) = log(P_P0(i))*(-tau(i));
end



%figure(2)
%%plot(t_isothermal,P_P0)
%title("$P_{tank}/P_{0}$ vs Time","Interpreter","latex");
%xlabel("Time [seconds]");
%ylabel("$P_{tank}/P_{0}$","Interpreter","latex");
%grid on 

%Finding A_star for current Supersoinc Tunnel 
A = 0.0232; %m^2
[M,T_ratio,P_ratio,rho_ratio,A_ratio] = flowisentropic(gamma,M);

A_star_tunnel = A/A_ratio;
tau_tunnel = (V_tank*((gamma+1)/2)^((gamma+1)/(2*(gamma-1))))/(C_d*A_star_tunnel*a_0);

t_isothermal_tunnel = [];
for i = 1:numel(P_P0)
    t_isothermal_tunnel(i) = log(P_P0(i))*(-tau_tunnel);
end

figure(3)
plot(t_isothermal_tunnel,P_P0);
title("$P_{tank}/P_{0}$ vs Time","Interpreter","latex");
xlabel("Time [seconds]","Interpreter","latex");
ylabel("$P_{tank}/P_{0}$","Interpreter","latex");
grid on 




%%Finding P_tank/p_0 for wind tunnel 
%Page 73 of High speed wind tunnel testing textbook 
%P_tank_tunnel = 448159; %initial tank pressure 
%P_tank_tot = (1/P_ratio)*P_tank_tunnel;%Total pressure of tank
%t_norm = t_tunnel*((A_star_tunnel*T_0*)
P_i = 689476; 
P_t = 448159;
P_f = P_t*1.5;
A_star_ft = A_star_tunnel*10.7639;
V_tank_ft = (V_tank*1550)/144;
T_0_far = (T_0 -273.15)*(9/5) + 32;
P_i_psi = P_i*0.000145038;
P_t_psi = P_t* 0.000145038;
P_f_psi = P_t_psi*1.5;
P_t_psf = P_t_psi *144;
T_0_rank = T_0_far + 460;
n = 1;
t_tunnel = 0.0353*(V_tank_ft/A_star_ft)*((sqrt(T_0_rank)/T_0_rank)*(P_i_psi/P_t_psi)*(1-(P_f_psi/P_i_psi)^(1/n)))




%% P_compressor vs time 

%Constants
M = 2.5;
P_0 = 448159; %pa
V_tank = 10.3; %m^3
gamma = 1.4;
C_d = 0.995;
A = 0.0232; %m^2
[M,T_ratio,P_ratio,rho_ratio,A_ratio] = flowisentropic(gamma,M);
A_star_tunnel = A/A_ratio;
T_0 = 319; %kelivn
R = 287;
a_0 = ((gamma*R*T_0)/M)^(1/2);

tau_tunnel = (V_tank*((gamma+1)/2)^((gamma+1)/(2*(gamma-1))))/(C_d*A_star_tunnel*a_0);
t_list = [0:1:15];
t_char_list = t_list./tau_tunnel;

P_tank = [];
for i = 1:numel(t_list)
    P_tank(i) = (P_0*(exp(-t_char_list(i))))/1000;
end

figure(4)
plot(t_list,P_tank)
grid on
xlabel("Time [seconds]","Interpreter","latex")
ylabel("$P \:[kPa]$","Interpreter","latex")
title("Tank Pressure vs Time","Interpreter","latex");

% Mass flow vs Time

rho_0 = P_0/(R*T_0);
rho_tank = [];
mass_flow = [];
for i = 1:numel(t_char_list)
    rho_tank(i) = rho_0*exp(-t_char_list(i));
end

for i = 1:numel(rho_tank)
    mass_flow(i) = C_d*A_star_tunnel*rho_tank(i)*sqrt(((gamma*R*T_0)/M))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    
end

figure(5)
plot(t_list,mass_flow)
grid on
xlabel("Time [seconds]","Interpreter","latex");
ylabel("$\dot{m}\: [Kg/s]$","Interpreter","latex");
title("Mass Flow vs Time","Interpreter","latex");

%Transient flow anlaysis assuming adiabitc condition 
rho_tank_ad = [];
P_tank_ad = [];
mass_flow_ad = [];

for i = 1:numel(t_char_list)
    rho_tank_ad(i) = rho_0*(1+((gamma-1)/2)*t_char_list(i))^(2/(1-gamma));
    P_tank_ad(i) = (P_0*(1+((gamma-1)/2)*t_char_list(i))^((2*gamma)/(1-gamma)))/1000;
    T_tank_ad(i) = T_0*(1+((gamma-1)/2)*t_char_list(i))^(-2);
    mass_flow_ad(i) = C_d*A_star_tunnel*rho_tank_ad(i)*sqrt(((gamma*R*T_tank_ad(i))/M))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
end
figure(6)
plot(t_list,P_tank_ad)
xlabel("T [sec]","interpreter","latex");
ylabel("$P_{tank} \: [kPa]$","interpreter","latex");
title("Tank Pressure vs Time [Adiabatic]","interpreter","latex");
grid on 

figure(7) 
plot(t_list,mass_flow_ad)
xlabel("T  [sec]","interpreter","latex");
ylabel("$\dot{m} \: [kg/s]$","interpreter","latex");
title("Mass flow vs Time [Adiabatic]","interpreter","latex");
grid on 

%% P tank calc quick
M = 2.5;
P_0 = 448159; %pa
V_tank = 10.3; %m^3
gamma = 1.4;
C_d = 0.995;
A = 0.0232; %m^2
[M,T_ratio,P_ratio,rho_ratio,A_ratio] = flowisentropic(gamma,M);
A_star_tunnel = A/A_ratio;
T_0 = 319; %kelivn
R = 287;
a_0 = ((gamma*R*T_0)/M)^(1/2);

tau_tunnel = (V_tank*((gamma+1)/2)^((gamma+1)/(2*(gamma-1))))/(C_d*A_star_tunnel*a_0);
t_list = [0:1:15];
t_char_list = t_list./tau_tunnel;
rho_0 = P_0/(R*T_0);
rho_tank = [];
mass_flow = [];
for i = 1:numel(t_char_list)
    rho_tank(i) = rho_0*exp(-t_char_list(i));
end

t_list = []

for i = 1:numel(rho_tank)
    mass_flow(i) = C_d*A_star_tunnel*rho_tank(i)*sqrt(((gamma*R*T_0)/M))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    t_list = (V_tank*rho_tank(i))/mass_flow(i)
end
