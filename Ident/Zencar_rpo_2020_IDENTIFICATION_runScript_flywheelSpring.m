clc
clear all
close all

load('Zencar_rpo_2020_IDENTIFICATION_spesession_flywheelSpring.mat')

cfi=SDOSessionData.Data.Workspace.LocalWorkspace.Final_estimation.Parameters(1).Value;
ratio=SDOSessionData.Data.Workspace.LocalWorkspace.Final_estimation.Parameters(2).Value;

U=12;
tsim=40;
tchange=12;

%starting parameters for estimation
% I=0.0467;
% k=0.312;
% b=0.03725;
% 
% R=17.1102;
% L=(6.0*10^-6);
% 
% cfi=0.1;
% ratio=0.1;

%final parameters with ratio and cfi loaded from spesession
I=0.0467*ratio;
k=0.312*ratio;
b=0.03725*ratio;

R=17.1102;
L=(6.0*10^-6);
cfi=cfi;


%freewheel model
A=[0 1 ; -k/I -b/I];
B=[0 ; 1/I];
C=[1 0];
D=[0];

%whole system model
Asys=[0 1 0 ; -k/I -b/I cfi/I ; 0 -cfi/L -R/L];
Bsys=[0 ; 0 ; 1/L];
Csys=[1 0 0;0 0 1];
Dsys=[0;0];

%model for parameter estimation
% Asys=[0                         1                          0;
%      -(k*m_ratio)/(I*m_ratio)   -(b*m_ratio)/(I*m_ratio)   cfi/(I*m_ratio);
%      0                          -cfi/L                     -R/L];
% Bsys=[0 ; 0 ; 1/L];
% Csys=[1 0 0;0 0 1];
% Dsys=[0;0];

out=sim('Zencar_rpo_2020_IDENTIFICATION_simModel_flywheelSpring.slx');

%% freewhell with nonzero initial angle
t=out.ficompare.time;
sysfi=out.ficompare.signals.values(:,1);
%shift signal
tcut=1632;
sysfi=sysfi(tcut:length(sysfi));
sysfi=[sysfi;zeros(length(t)-length(sysfi),1)];

estfi=out.ficompare.signals.values(:,2);

tend=12;

figure
plot(t(1:min(find(t>=tend))),sysfi(1:min(find(t>=tend))),'k','linewidth',2);hold on
plot(t(1:min(find(t>=tend))),estfi(1:min(find(t>=tend))),'--r','linewidth',2)
grid on
title('Pure mechanical model , parameter estimation')
legend('measured angle','simulated angle')
xlabel('time [s]')
ylabel('system angle[rad]')
xlim([0,tend])
set(gcf,'units','normalized','outerposition',[0,0,1,1]);

%% Whole system
t=out.syscompare.time;
sysfi=out.syscompare.signals.values(:,1);
syscur=out.syscompare.signals.values(:,3);
estfi=out.syscompare.signals.values(:,2);
estcur=out.syscompare.signals.values(:,4);

tend=25;

figure
subplot(2,1,1)
plot(t(1:min(find(t>=tend))),sysfi(1:min(find(t>=tend))),'-k','linewidth',2);hold on
plot(t(1:min(find(t>=tend))),estfi(1:min(find(t>=tend))),'--r','linewidth',2)
grid on
title('Electromechanical model , parameter estimation')
legend('measured angle','simulated angle')
xlabel('time [s]')
ylabel('system angle[rad]')
xlim([0,tend])

subplot(2,1,2)
lim=length(t);
plot(t(1:min(find(t>=tend))),syscur(1:min(find(t>=tend))),'-k','linewidth',2);hold on
plot(t(1:min(find(t>=tend))),estcur(1:min(find(t>=tend))),'--r','linewidth',2)
grid on
legend('measured current','simulated current')
xlabel('time [s]')
ylabel('motor current[A]')
xlim([0,tend])

set(gcf,'unit','normalized','outerposition',[0,0,1,1])
