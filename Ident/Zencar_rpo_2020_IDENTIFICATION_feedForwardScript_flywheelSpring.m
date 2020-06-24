clc
clear all
close all

load('Zencar_rpo_2020_IDENTIFICATION_spesession_flywheelSpring.mat')

cfi=SDOSessionData.Data.Workspace.LocalWorkspace.Final_estimation.Parameters(1).Value;
ratio=SDOSessionData.Data.Workspace.LocalWorkspace.Final_estimation.Parameters(2).Value;

%final parameters with ratio and cfi loaded from spesession
I=0.0467*ratio;
k=0.312*ratio;
b=0.03725*ratio;

R=17.1102;
L=(6.0*10^-6);
cfi=cfi;

%whole system model
Asys=[0 1 0 ; -k/I -b/I cfi/I ; 0 -cfi/L -R/L];
Bsys=[0 ; 0 ; 1/L];
Csys=[1 0 0;0 0 1];
Dsys=[0;0];

fi_wanted=2;
out=sim('Zencar_rpo_2020_IDENTIFICATION_feedForwardModel_flywheelSpring.slx')

%% plotting
t=out.feedForward.time;
fi=out.feedForward.signals.values(:,1);

figure
plot(t,fi,'b','linewidth',2)
grid minor
title('FeedForward regulation')
legend('angle')
xlabel('time [s]')
ylabel('position [rad]');

set(gcf,'units','normalized','outerposition',[0 0 1 1])