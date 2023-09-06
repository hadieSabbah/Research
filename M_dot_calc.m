%% Calculating mass flow 
P0 = [450:20:690]; %kpa
T0 = 300; %kelvin %Assuming constant stagnation temperature 
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
xlabel("Stagnation Pressure[kPa]")
ylabel("Mass Flow[Kg/s]")
title("Mass Flow Vs Stagnation Pressure")
grid on






