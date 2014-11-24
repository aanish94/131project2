clc
clear all
close all

%Defining Problem Parameters

%(meters)
D = 2.4384;
R = D / 2;
L = 9.144;
t_al = 5/1000;
% t_cork; %we are solving for this

%if ho = 25, === t_cork = 0.35 mm
%if ho = 5 === t_cork = 8 mm

%(W m^-2 K^-1)
ho = 25;
hi = 50;

%(Wm^-1K^-1
k_cork = 0.07;
k_al = 205;

%in Kelvin
Tinf = 300.15;
Tox = 90;
Tso = 175;

sigma = 5.67*10^-8;
epsilon = 0.7;
% A = pi * D * L;
% q_rad = sigma * epsilon * A * ( Tinf^4 - Tso^4);
%
% q_conv = ho * A * (Tinf - Tso);
%
% q_total = q_rad + q_conv
%Liquid Oxygen enthalpy of formation @ 90 K  ( kJ/ kg)
hfg =  212.3;
%% Calculating Thermal resistances
figure(1)
for ho = 5:1:25
    difference = 100;
    t_cork_final = 0;
    for t_cork = 0:0.0001/1000:10/1000 %0.353/1000
        %Outer Convection Resistance
        D_outer = 2 * (R + t_al + t_cork);
        Ao = pi * D_outer * L;
        
        Rconv_1 = 1 / (ho*Ao);
        
        %Cork Conduction Resistance
        r_outerc = R + t_al + t_cork;
        r_innerc = R + t_al;
        
        Rcond_cork = log(r_outerc / r_innerc) / (2*pi*L*k_cork);
        
        %Aluminum Conduction Resistance
        r_outeral = R + t_al;
        r_inneral = R;
        
        Rcond_al = log(r_outeral / r_inneral) / (2*pi*L*k_al);
        
        %Inner Convection Resistance
        D_inner = D;
        Ai = pi * D_inner * L;
        
        Rconv_2 = 1 / (hi * Ai);
        
        %Solving Equation
        RHS = (Tso - 90) / ( Rcond_cork + Rcond_al + Rconv_2);
        
        q_rad = sigma * epsilon * Ao * (Tinf^4 - Tso^4);
        q_conv = (Tinf - Tso) / Rconv_1;
        LHS = q_rad + q_conv;
        if abs(LHS - RHS) < difference
            t_cork_final = t_cork;
            difference = abs(LHS - RHS);
        end
    end
    q_total = q_rad + q_conv;
    m_dot = q_total / ( hfg * 10^3);
    
    %% PART B
    figure(1)
    plot(ho,t_cork_final*1000,'x')
    hold on
    xlabel('Outside heat transfer coefficient (W/m^2K)','Fontsize',14)
    ylabel('Cork Thickness (mm)','Fontsize',14)
    title('Effect of outside heat transfer on cork thickness','Fontsize',16)
    %Convert to mm
    fprintf('Thickness of cork: %5.2f mm\n',t_cork_final*1000)
    
    %% PART A
    figure(2)
    plot(ho,q_total,'x')
    hold on
    xlabel('Outside heat transfer coefficient (W/m^2K)','Fontsize',14)
    ylabel('q total (W)','Fontsize',14)
    title('Effect of outside heat transfer on total heat transfer','Fontsize',16)
    fprintf('Heat transfer to LOX, qtotal = %5.2f W\n',q_total)
    
    %% PART C
    fprintf('Mass flow rate, m_dot = %5.2f kg/s\n',m_dot)
    figure(3)
    plot(ho,m_dot,'x')
    hold on
    xlabel('Outside heat transfer coefficient (W/m^2K)','Fontsize',14)
    ylabel('m_dot (kg/s)','Fontsize',14)
    title('Effect of outside heat transfer coefficient on mass flow rate','Fontsize',16)
    
    disp(' ' )
end
