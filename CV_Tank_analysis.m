%Problem Analysis 
clear; clc;
R = 287; %j/kg*k
gamma = 1.4;

P_i = 101325; %Pa
P_tank = [];
P_f = 1.31e+6; %Pa, The final pressure of the tank 
T_i = 300; %Kelvin
T_f = T_i; %Assuming Isothermal. 
r = 1/2; %m. This is the radius of the tank
h = 4; %m. This is the height of the tank
V_tank = pi*h*r^2; %Volume of Tank

rho_i = P_i/(R*T_i); %kg/m^3
rho_f = P_f/(R*T_f); %kg/m^3
rho_ambient = 1.225; %kg/m^3

%%%% Transient Analysis based on continuity Equation %%%%

t_list = [];
rho_list = linspace(1,rho_f/rho_ambient,30); %list used for flowisentropic tables

m_in = 3; %kg/s. This is the mass flow rate going into the tank
A_star_list = linspace(0.0088,1,numel(rho_list));
A_star = 0.0088*5;
%Initiliazing Variables
rho_tank = [];
P_tank = [];
m_out = [];
M = [];
%Finding the time step for each density ratio

%%% CREATE A CELL ARRAY FOR EACH A_STAR_LIST VALUE TO DETERMINE THE
%%% PROPERTIES OF THESE VALUES AT DIFFERENT AREAS. THIS CRITICAL FOR THE
%%% SIZING STAGE OF THE TUNNEL. 
%for n = numel(A_star_list)
    for i = 1:numel(rho_list)
        [M(i),T,P,rho_ratio,area] = flowisentropic(gamma,(1/rho_list(i)),'dens');
        rho_tank(i) = (1/rho_ratio)*rho_ambient;
        P_tank(i) = rho_tank(i)*R*T_i;
            if M(i) > 1
                M(i) = 1;
                m_out(i) = (((A_star*P_tank(i)))/(T_i))*sqrt(gamma/R)*(((gamma+1)/2))^((-gamma-1)/(2*(gamma-1)));
            else
                m_out(i) = (((A_star*P_tank(i)))/(T_i))*sqrt(gamma/R)*M(i)*(1 + ((gamma-1)/2)*M(i)^2)^((-gamma-1)/(2*(gamma-1)));
            end
    
        t_list(i) = (m_out(i) - rho_i*V_tank +rho_tank(i)*V_tank)/(m_in);

    
    end
%end

%Figures

%Mass flow out vs Time
figure(1)
plot(t_list,m_out)
xlabel("Time")
ylabel("Exit Mass flow[Kg/s]")
title("Exit Mass Flow Vs Time")
grid on

%Pressure in the tank vs Time
figure(2)
plot(t_list,P_tank./1000)
xlabel("Time")
ylabel("$P_{tank} [kPa]$","Interpreter","latex");
title("Tank Pressure Vs Time")
grid on


figure(3)
plot(t_list,M)
xlabel("Time[S]")
ylabel("Mach")
title("Mach vs Time")
grid on 

%% Inflow Tank Analysis %% 
t_inflow = [];
for i = 1:numel(rho_tank)

    t_inflow(i) = (V_tank*rho_tank(i) - V_tank*rho_i)/(m_in);
    
end

figure(4)
plot(t_inflow,P_tank/1000)
xlabel("Time[S]")
ylabel("Tank Pressure [kPa]")
title("Tank Pressure Vs Time")
grid on 
hold on 
plot(t_list,P_tank./1000)
xlabel("Time")
ylabel("$P_{tank} [kPa]$","Interpreter","latex");
title("Tank Pressure Vs Time")
grid on
legend("Closed Vessel","Open Vessel","location","best")


