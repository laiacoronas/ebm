%% OXYGEN DIFFUSION
% Marta Maldonado and Cristina Mendoza

%% OXYGEN DIFFUSION - MEMBRANE 
% 1. Simulate the different experiments in table 1. (Vr = 50ml, C0donnor = 100%, C0receptor = 0%)

t = linspace(0,2e5,1000);

%%
%SIS-L:
N_sisl = 14;
z_sisl = 48.64;
b_sisl = 13.39;
D_sisl = 2.43e-6;

C_sisl = concentration_oxygen(100,D_sisl,t,b_sisl);

figure(1)
plot(t,C_sisl)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%%
%SIS_A:
N_sisa = 12;
z_sisa = 48.64;
b_sisa = 13.39;
D_sisa = 5.63e-6;

C_sisa = concentration_oxygen(100,D_sisa,t,b_sisa);

figure(2)
plot(t,C_sisa)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%%
%UBM-L:
N_ubml = 9;
z_ubml = 99.54;
b_ubml = 5.69;
D_ubml = 4.64e-6;

C_ubml = concentration_oxygen(100,D_ubml,t,b_ubml);

figure(3)
plot(t,C_ubml)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%%
%UBM-A:
N_ubma = 12;
z_ubma = 99.54;
b_ubma = 5.69;
D_ubma = 3.75e-6;

C_ubma = concentration_oxygen(100,D_ubma,t,b_ubma);

figure(4)
plot(t,C_ubma)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%%
%UBS:
N_ubs = 16;
z_ubs = 111.81;
b_ubs = 6.48;
D_ubs = 6.61e-6;

C_ubs = concentration_oxygen(100,D_ubs,t,b_ubs);

figure(5)
plot(t,C_ubs)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%%
%Dacron:
N_dacron = 12;
z_dacron = 99.54;
b_dacron = 5.69;
D_dacron = 3.75e-6;

C_dacron = concentration_oxygen(100,D_dacron,t,b_dacron);

figure(6)
plot(t,C_dacron)
title('Concentration of dissolved oxygen in the reciever chamber')
ylabel('CR(t) [mg/l]')
xlabel('Time [s]')


%% 2. Analyze how the different parameters of the membrane impacts the effective time constant

% If z (thickness) increases, the time constant (tau) also increases. 
%However, if the diffusion coefficient increases, the time constant
%decreases.
%It can also be observed that the area increases as the time constant
%decreases.

Vr = 50;
A_sisl = linspace(1,5,100);
T_sisl = 1./(A_sisl*D_sisl/z_sisl/Vr)*10e-9;

figure(7)
plot(A_sisl,T_sisl);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')

%%
A_sisa = linspace(1,5,100);
T_sisa = 1./(A_sisa*D_sisa/z_sisa/Vr)*10e-9;
figure(8)
plot(A_sisa,T_sisa);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')

%%
A_ubml = linspace(1,5,100);
T_ubml = 1./(A_ubml*D_ubml/z_ubml/Vr)*10e-9;
figure(9)
plot(A_ubml,T_ubml);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')


%%
A_ubma = linspace(1,5,100);
T_ubma = 1./(A_ubma*D_ubma/z_ubma/Vr)*10e-9;
figure(10)
plot(A_ubma,T_ubma);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')


%%
A_ubs = linspace(1,5,100);
T_ubs = 1./(A_ubs*D_ubs/z_ubs/Vr)*10e-9;
figure(11)
plot(A_ubs,T_ubs);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')


%%
A_dacron = linspace(1,5,100);
T_dacron = 1./(A_dacron*D_dacron/z_dacron/Vr)*10e-9;
figure(12)
plot(A_dacron,T_dacron);
title('Time constant vs Diffusion area')
xlabel('Diffusion area of membrane [cm^{2}]')
ylabel('Time constant [s^{-1}] *10^{4}')

%% OXYGEN DIFFUSION - CYLINDER
%1. Determine C(t) for four different cylinder diameters: 0.5, 1, 2 and 4 cm for a given change
%of oxygen concentration and D fixed at 4.5x10-6 cm2/s. Estimate the "effective” time
%constant (time to reach a certain percentage of the total value of oxygen concentration). Is
%the "effective” time constant linear with the diameter?


diameter = [0.5,1,2,4];   %cm
r = diameter./2;  %cm
D = 4.5e-6; %cm2/s
t = linspace(0,3e5,100);   %s

%first 3 roots of the Bessel functions:
alpha_1 = 2.4048;
alpha_2 = 5.5201;
alpha_3 = 8.6537;

%Bessel functions:
J_1 = besselj(1,alpha_1);
J_2 = besselj(1,alpha_2);
J_3 = besselj(1,alpha_3);
  
%Concentration C(t)
for i = 1:4
    for a = 1:length(t)
        T(i,a) = D.*t(a)/r(i)^2;

        alpha_list = [alpha_1,alpha_2,alpha_3];
        
        suma = 0;

        if T(i,a) > 0.06301
            
            for n = 1:3
                f = (exp(-alpha_list(n)^2*T(i,a)))/(alpha_list(n)*(besselj(1,alpha_list(n))));
                suma = f+suma;
            end
            C(i,a) = 1-2.*suma;
        

        elseif T(i,a) == 0.06301
           C(i,a) = 0.03575;
        

        elseif T(i,a) < 0.06301
            %Fractional Bessel functions:
            N(i,a) = 1./(8.*T(i,a));
            K_1_4(i,a) = besselk(0.25,N(i,a));
            K_3_4(i,a) = besselk(0.75,N(i,a));

            C(i,a) = (exp(-1./(8.*T(i,a)))./sqrt(pi.*T(i,a))).*(K_1_4(i,a).*(1./(8.*T(i,a)))-T(i,a)./4.*K_3_4(i,a).*(1./(8.*T(i,a))));

            %C(i,a) = (exp(-1./(8.*T(i,a)))./sqrt(pi.*T(i,a))).*(K_1_4(i,a))-(T(i,a)./4).*K_3_4(i,a);
        end 
    end
   
end
C1 = C*160;  %mmHg

plot(t(1:length(C1)),C1)
ylabel("Oxygen partial pressure (mmHg)")
xlabel("Time (s)")
title("Cylinder diameter effect on oxygen diffusion")
%yline(160*0.5)
legend(["0.5 cm","1 cm", "2 cm", "4 cm"])

%% Estimate the effective time constant (time to reach a certain percentage
% of the total value of oxygen concentration). Is the "effective" time
% constant linear with the diameter?

%T = D*/t*r^2.
% It is not linear with the diameter because it depends on the inverse of
% the power of 2 of the diameter (or radius).

perc = 0.2;   %we choose to study the time constant at oxygen concentration of 50%
C_end = 160;

C_perc = C_end*perc;

for i = length(C1(1,:))
    if C1(:,i) == C_perc
        T_const = T(i,:);
    end
end

%no hem aconseguit res amb el codi anterior

%% 2. Same as task 1, but fixing te diameter (2 cm) and testing the different difussion coefficients in table 1. Is the effective time constant linear with the difussion coefficient?
d = 2;
rad = d/2;
D_coef = [4.5e-6, 3.5e-6, 1.7e-6];

for i = 1:3
    T2(i,:) = D_coef(i).*t/rad^2;

    %first 3 roots of the Bessel functions:
    alpha_1 = 2.4048;
    alpha_2 = 5.5201;
    alpha_3 = 8.6537;
    
    %Bessel functions:
    J_1 = besselj(1,alpha_1);
    J_2 = besselj(1,alpha_2);
    J_3 = besselj(1,alpha_3);
    
  
    %Concentration C(t)
    for a = 1:length(t)

        %Fractional Bessel functions:
        N2(i,a) = 1./(8.*T2(i,a));
        K_1_4_2(i,a) = besselk(0.25,N2(i,a));
        K_3_4_2(i,a) = besselk(0.75,N2(i,a));
        
        alpha_list = [alpha_1,alpha_2,alpha_3];

        if T2(i,a) < 0.06301
            C2(i,a) = (exp(-1./(8.*T2(i,a)))./sqrt(pi.*T2(i,a))).*(K_1_4_2(i,a).*(1./(8.*T2(i,a)))-T2(i,a)./4.*K_3_4_2(i,a).*(1./(8.*T2(i,a))));
        end

        if T2(i,a) > 0.06301
            suma = 0;
            for n = 1:3
                f2 = (exp(-alpha_list(n)^2*T2(i,a)))/(alpha_list(n)*(besselj(1,alpha_list(n))));
                suma = f2+suma;
            end
            C2(i,a) = 1-2.*suma;
        end

        if T2(i,a) == 0.06301
           C2(i,a) = 0.03575;
        end
    end
end
C3 = C2*160;  %mmHg

plot(t(1:length(C3)),C3)
ylabel("Oxygen partial pressure (mmHg)")
xlabel("Time (s)")
title("Different diffusion coefficients effect on oxygen diffusion")
legend(["4.5e-6","3.5e-6","1.7e-6"])


% The effective time constant is linear with the diffusion coefficient
% because T = D*t/r^2. Therefore there is a linear relationship between T
% and D.

%% Optional: Fit the data simulated to an exponential fitting. Is the fitting good? In which range?
