%% Oxygen diffusion membrane

%% Task 1
z= [48.64, 48.64, 99.54, 99.54, 111.81, 630.83];
beta= [13.39, 13.39, 5.69, 5.69, 6.48, 0.95];
dprima=[2.43, 5.63, 4.64, 3.75, 6.61, 46.67];
t=1:1000:2*10^5;
Cr=[];
a=[];
for i=1:6
    for j=1:200
        a(i,j)=exp(-beta(i)*t(j)*dprima(i)*10^(-6));
        Cr(i, j)=(100*(1-a(i,j)))/(2-a(i,j));
    end
    plot(t,Cr(i,:))
    hold on
end
hold off

%% Task 2

A=1:0.01:5;
tau=[];
Vr=50
for i=1:6
    for j=1:401
        tau(i,j)=1/(A(j)*z(i)*10^(-4)*(1/Vr)*dprima(i)*10^(-6));
    end
    plot(A, tau(i,:))
    hold on
end
hold off

%%
d = [0.5, 1, 2, 4];
D = 4.5*10^(-6);
t = 0:1:3*10^5;
r = d./2;
T = [];
for i = 1:length(r)
    for j = 1:length(t)
        T(i, j) = D*t(j)/r(i)^2;
    end
    plot(t, T(i, :))
    hold on
end
hold off

% It is apparently linnear!!

%%

alpha_1 = 2.4048; 
alpha_2 = 5.5201;
alpha_3 = 8.6537;
alpha = [2.4048, 5.5201, 8.6537];

J_1 = besselj(1,alpha_1);
J_2 = besselj(1,alpha_2);
J_3 = besselj(1,alpha_3);
J = [J_1, J_2, J_3];

%N = 1./(8*T);
%K_1_4 = besselk(0.25,N);
%K_3_4 = besselk(0.75,N);

sumatori = 0;
C = [];
threshold = 0.06301;

for i = 1:length(r)
    for j = 1:length(t)
        T = (D*t(j))/r(i)^2;
        sumatori = 0;
        if T > threshold
            for m = 1:3
                sumatori = sumatori + exp(-(alpha(m)^2)*T)/(alpha(m)*J(m));
            end
            C(i,j) = 1 - 2*sumatori;
        elseif T == threshold
            C(i, j) = 0.03575;

        elseif T < threshold
            N = 1/(8*T);
            K_1_4 = besselk(0.25,N);
            K_3_4 = besselk(0.75,N);
            C(i, j) = (exp(-N)/sqrt(pi*T))*(K_1_4-(T/4)*K_3_4);
        end
    end
    plot(t, C(i,:))
    hold on
end
xlabel("Time [seconds]")
ylabel("Oxygen Partial Pressure")
title("Cylinder Diameter effect on Oxygen Diffusion")
legend('0,5 cm', '1 cm', '2 cm', '4 cm')
hold off

%%

D = [4.5 * 10^-6, 3.5*10^-6, 1.7*10^-6]
d = 2
r = d/2;

alpha_1 = 2.4048; 
alpha_2 = 5.5201;
alpha_3 = 8.6537;
alpha = [2.4048, 5.5201, 8.6537];

J_1 = besselj(1,alpha_1);
J_2 = besselj(1,alpha_2);
J_3 = besselj(1,alpha_3);
J = [J_1, J_2, J_3];

%N = 1./(8*T);
%K_1_4 = besselk(0.25,N);
%K_3_4 = besselk(0.75,N);

sumatori = 0;
C = [];
threshold = 0.06301;

for i = 1:length(D)
    for j = 1:length(t)
        T = (D(i)*t(j))/r^2;
        sumatori = 0;
        if T > threshold
            for m = 1:3
                sumatori = sumatori + exp(-(alpha(m)^2)*T)/(alpha(m)*J(m));
            end
            C(i,j) = 1 - 2*sumatori;
        elseif T == threshold
            C(i, j) = 0.03575;

        elseif T < threshold
            N = 1/(8*T);
            K_1_4 = besselk(0.25,N);
            K_3_4 = besselk(0.75,N);
            C(i, j) = (exp(-N)/sqrt(pi*T))*(K_1_4-(T/4)*K_3_4);
        end
    end
    plot(t, C(i,:))
    hold on
end
xlabel("Time [seconds]")
ylabel("Oxygen Partial Pressure")
title("Cylinder Diameter effect on Oxygen Diffusion")
legend('4,5*10^6 cm', '3,5*10^6 cm', '1,7*10^6 cm')
hold off


%%

rng('default') ;
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
x = 1:1:length(t);
y = stress(:,i);
b0 = [0,1]; % initial parameters for the routine to start the fitting
modelstr_prova = 'y ~ b1*exp(b2*x)-b1';
mdl_prova = fitnlm(x,y,modelstr_prova,b0,'Options',opts);
alpha(1,i)=table2array(mdl_prova.Coefficients(2,1));
beta(1,i)=table2array(mdl_prova.Coefficients(1,1));



