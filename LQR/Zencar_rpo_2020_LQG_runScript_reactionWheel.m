clc
clear all
close all
warning('off')

%% parameters and settings init

simulate=0;% simulink simulation or load from memory , if 0 overwrites all settings

visualize=1;% if 1 System will be simulated in figure (make sure plotLinear is set right)
plotLinear=0;% 0 if using nonlinear plant , 1 if using linear plant


fi0=-0.20;%initial pendulum angle


g=9.81;

bw=0.1;%viscous friction of pendulum arm joint
bm=0.05;%viscous friction of motor

rm=0.055;%half of motor diameter
mm=0.1;%mass of motor
lm=0.4;%distance from 0 to motor

rw=0.2;%half of wheel diameter
mw=0.4;%wheel mass
lw=0.4;% distance from 0 to wheel center


mb=0.1;%pendulum arm mass
lb=0.2;%distance from 0 to center of gravity of pendulum wheel arm;
rb=2*lb;%length of pendulum arm

Kmgl=(mm*lm + mw*mw + mb*lb)*g;%gravitational constant

Im=0.5*mm*rm^2;%moment of inertia of motor
Imo=Im+mm*lm^2;%moment of inertia of motor in respect to 0

Iw=mw*rw^2;%moment of inertia of wheel
Iwo=Iw+mw*lw^2;%moment of inertia of whell in respect to 0

Ibo=1/3*mb*rb^2;%moment of inertia of pendulum arm

I=Imo+Iwo+Ibo;%total moment of inertia

% zero degrees is up , linearization with sin(x)=x , positive angle is to the right
A= [0                     1                 0;
    +Kmgl/I               -bw/I             +bm/I; 
    -Kmgl/I               +bw/I             -(Iw+I)/(Iw*I)*bm];

Aadd=[0;0;0];

B= [0;
    -1/I;
    (Iw+I)/(Iw*I)];

C = [1  0  0;0 1 0; 0 0 1];
     
D=[0;0;0];





%% LQR init
Q=diag([50,10,5e-4]);
R=1;
K=lqr(A,B,Q,R);



%% Kalman init

sysc=ss(A,B,C,D)
ts=1e-2;
sysd=c2d(sysc,ts);
Ad=sysd.A;
Bd=sysd.B;

H=[1 0 0];
V = 1*eye(3);
W=  100;
xhat=[fi0+(-1)^randi(10)*rand*2*pi/3200,0,0]'; % initial guess is current angle plus noise
P=zeros(3);

disp('Linearized system poles: ')
poles=pole(sysc)

if(simulate)
    out=sim('Zencar_rpo_2020_LQG_simModel_reactionWheel.slx');
else
    load('Zencar_rpo_2020_LQG_simulationdata_reactionWheel.mat')
end
axis_headroom=0.1;

%% Open loop response plotting 
t=out.open_loop.time;
fi=out.open_loop.signals.values(:,1);
fidot=out.open_loop.signals.values(:,2);
omega=out.open_loop.signals.values(:,3);

figure
subplot(3,1,1)
plot(t,fi,'-b','linewidth',2)
grid on
legend('pendulum angle')
xlabel('time [s]')
ylabel('pendulum angle [rad]')
title('Open loop response')

subplot(3,1,2)
plot(t,fidot,'-b','linewidth',2)
grid on
legend('pendulum omega')
xlabel('time [s]')
ylabel('pendulum omega [rads-1]')

subplot(3,1,3)
plot(t,omega,'-b','linewidth',2)
grid on
legend('reaction wheel omega')
xlabel('time [s]')
ylabel('reaction wheel omega [rads-1]')


set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%% LQR zero setpoint with no noise plotting , no saturation
t=out.zero_reg.time;
fi=out.zero_reg.signals.values(:,1);
fidot=out.zero_reg.signals.values(:,2);
omega=out.zero_reg.signals.values(:,3);
action_value=out.zero_reg.signals.values(:,5);

tend=10;

t=t(1:tend/(t(length(t)))*length(t)+1);

figure
subplot(4,1,1)
plot(t,fi(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(fi)-(max(fi)-min(fi))*axis_headroom,max(fi)+(max(fi)-min(fi))*axis_headroom])
legend('pendulum angle')
xlabel('time [s]')
ylabel('pendulum angle [rad]')
title('LQR regulation to zero, weight on pendulum angle')

subplot(4,1,2)
plot(t,fidot(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(fidot)-(max(fidot)-min(fidot))*axis_headroom,max(fidot)+(max(fidot)-min(fidot))*axis_headroom])
legend('pendulum omega')
xlabel('time [s]')
ylabel('pendulum omega [rads-1]')

subplot(4,1,3)
plot(t,omega(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(omega)-(max(omega)-min(omega))*axis_headroom,max(omega)+(max(omega)-min(omega))*axis_headroom])
legend('reaction wheel omega')
xlabel('time [s]')
ylabel('reaction wheel omega [rads-1]')


subplot(4,1,4)
plot(t,action_value(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(action_value)-(max(action_value)-min(action_value))*axis_headroom,max(action_value)+(max(action_value)-min(action_value))*axis_headroom])
legend('action value')
xlabel('time [s]')
ylabel('Torque (Nm)')


set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%% LQR zero setpoint with no noise plotting ,saturation
t=out.zero_reg_sat.time;
fi=out.zero_reg_sat.signals.values(:,1);
fidot=out.zero_reg_sat.signals.values(:,2);
omega=out.zero_reg_sat.signals.values(:,3);
action_value=out.zero_reg_sat.signals.values(:,5);

tend=10;

t=t(1:tend/(t(length(t)))*length(t)+1);

figure
subplot(4,1,1)
plot(t,fi(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(fi)-(max(fi)-min(fi))*axis_headroom,max(fi)+(max(fi)-min(fi))*axis_headroom])
legend('pendulum angle')
xlabel('time [s]')
ylabel('pendulum angle [rad]')
title('LQR regulation to zero, weight on pendulum angle and pendulum ang. velocity')

subplot(4,1,2)
plot(t,fidot(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(fidot)-(max(fidot)-min(fidot))*axis_headroom,max(fidot)+(max(fidot)-min(fidot))*axis_headroom])
legend('pendulum omega')
xlabel('time [s]')
ylabel('pendulum omega [rads-1]')

subplot(4,1,3)
plot(t,omega(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(omega)-(max(omega)-min(omega))*axis_headroom,max(omega)+(max(omega)-min(omega))*axis_headroom])
legend('reaction wheel omega')
xlabel('time [s]')
ylabel('reaction wheel omega [rads-1]')


subplot(4,1,4)
plot(t,action_value(1:length(t)),'-b','linewidth',2)
grid on
axis([0,t(length(t)),min(action_value)-(max(action_value)-min(action_value))*axis_headroom,max(action_value)+(max(action_value)-min(action_value))*axis_headroom])
legend('action value')
xlabel('time [s]')
ylabel('Torque (Nm)')


set(gcf,'units','normalized','outerposition',[0 0 1 1]);


%% kalman plotting 
t=out.kalman_estimated.time;
fi_estimated=out.kalman_estimated.signals.values(:,1);
fidot_estimated=out.kalman_estimated.signals.values(:,2);
omega_estimated=out.kalman_estimated.signals.values(:,3);

fi_measured=out.kalman_measured.signals.values(:,1);

fi_real=out.kalman_real.signals.values(:,1);
fidot_real=out.kalman_real.signals.values(:,2);
omega_real=out.kalman_real.signals.values(:,3);



tend=6;

t=t(1:tend/(t(length(t)))*length(t)+1);

figure
subplot(3,1,1)
plot(t,fi_measured(1:length(t)),'-b','linewidth',1);hold on
plot(t,fi_estimated(1:length(t)),'-r','linewidth',2);
plot(t,fi_real(1:length(t)),'--g','linewidth',2);
grid on
axis([0,t(length(t)),min(fi_measured)-(max(fi_measured)-min(fi_measured))*axis_headroom,max(fi_measured)+(max(fi_measured)-min(fi_measured))*axis_headroom])
legend('measured angle','real angle','estimated angle')
xlabel('time [s]')
ylabel('pendulum angle [rad]')
title('Kalman filter estimation , W=50')

subplot(3,1,2)
plot(t,fidot_estimated(1:length(t)),'-r','linewidth',2);hold on
plot(t,fidot_real(1:length(t)),'--k','linewidth',2);
grid on
axis([0,t(length(t)),min(fidot_estimated)-(max(fidot_estimated)-min(fidot_estimated))*axis_headroom,max(fidot_estimated)+(max(fidot_estimated)-min(fidot_estimated))*axis_headroom])
legend('real pendulum velocity','estimated pendulum velocity')
xlabel('time [s]')
ylabel('pendulum omega [rads-1]')

subplot(3,1,3)
plot(t,omega_estimated(1:length(t)),'-r','linewidth',2);hold on
plot(t,omega_real(1:length(t)),'--k','linewidth',2);
grid on
axis([0,t(length(t)),min(omega_estimated)-(max(omega_estimated)-min(omega_estimated))*axis_headroom,max(omega_estimated)+(max(omega_estimated)-min(omega_estimated))*axis_headroom])
legend('real reaction wheel omega','estimated reaction wheel omega')
xlabel('time [s]')
ylabel('reaction wheel omega [rads-1]')


set(gcf,'units','normalized','outerposition',[0 0 1 1]);


%% simulation , nonlinear plant
if(plotLinear==0)
	t=out.nonlinear.time;
	fi=out.nonlinear.signals.values(:,1);
	fidot=out.nonlinear.signals.values(:,2);
	omega=out.nonlinear.signals.values(:,3);
	action_value=out.nonlinear.signals.values(:,5);

	tend=12;

	t=t(1:tend/(t(length(t)))*length(t)+1);

	figure
	subplot(4,1,1)
	plot(t,fi(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(fi)-(max(fi)-min(fi))*axis_headroom,max(fi)+(max(fi)-min(fi))*axis_headroom])
	legend('pendulum angle')
	xlabel('time [s]')
	ylabel('pendulum angle [rad]')
	title('System simulation with 2 external force impules (nonlinear plant)')

	subplot(4,1,2)
	plot(t,fidot(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(fidot)-(max(fidot)-min(fidot))*axis_headroom,max(fidot)+(max(fidot)-min(fidot))*axis_headroom])
	legend('pendulum omega')
	xlabel('time [s]')
	ylabel('pendulum omega [rads-1]')

	subplot(4,1,3)
	plot(t,omega(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(omega)-(max(omega)-min(omega))*axis_headroom,max(omega)+(max(omega)-min(omega))*axis_headroom])
	legend('reaction wheel omega')
	xlabel('time [s]')
	ylabel('reaction wheel omega [rads-1]')


	subplot(4,1,4)
	plot(t,action_value(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(action_value)-(max(action_value)-min(action_value))*axis_headroom,max(action_value)+(max(action_value)-min(action_value))*axis_headroom])
	legend('action value')
	xlabel('time [s]')
	ylabel('Torque (Nm)')


    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

end

%% simulation , linear plant
if(plotLinear==1)
	t=out.linear.time;
	fi=out.linear.signals.values(:,1);
	fidot=out.linear.signals.values(:,2);
	omega=out.linear.signals.values(:,3);
	action_value=out.linear.signals.values(:,5);

	tend=12;

	t=t(1:tend/(t(length(t)))*length(t)+1);

	figure
	subplot(4,1,1)
	plot(t,fi(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(fi)-(max(fi)-min(fi))*axis_headroom,max(fi)+(max(fi)-min(fi))*axis_headroom])
	legend('pendulum angle')
	xlabel('time [s]')
	ylabel('pendulum angle [rad]')
	title('System simulation with 2 external force impules (linear plant)')

	subplot(4,1,2)
	plot(t,fidot(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(fidot)-(max(fidot)-min(fidot))*axis_headroom,max(fidot)+(max(fidot)-min(fidot))*axis_headroom])
	legend('pendulum omega')
	xlabel('time [s]')
	ylabel('pendulum omega [rads-1]')

	subplot(4,1,3)
	plot(t,omega(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(omega)-(max(omega)-min(omega))*axis_headroom,max(omega)+(max(omega)-min(omega))*axis_headroom])
	legend('reaction wheel omega')
	xlabel('time [s]')
	ylabel('reaction wheel omega [rads-1]')


	subplot(4,1,4)
	plot(t,action_value(1:length(t)),'-b','linewidth',2)
	grid on
	axis([0,t(length(t)),min(action_value)-(max(action_value)-min(action_value))*axis_headroom,max(action_value)+(max(action_value)-min(action_value))*axis_headroom])
	legend('action value')
	xlabel('time [s]')
	ylabel('Torque (Nm)')


    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

end


%% visualization
if(visualize)
    figCur=figure;
    %figure(figCur);
    plt=plot(0,0);
    axis equal
    grid on
    if plotLinear
        time=out.linear.time;
        fi=out.linear.signals.values(:,1);
        fidot=out.linear.signals.values(:,2);
        omegaw=out.linear.signals.values(:,3);
        anglew=out.linear.signals.values(:,4);
    else
        time=out.nonlinear.time;
        fi=out.nonlinear.signals.values(:,1);
        fidot=out.nonlinear.signals.values(:,2);
        omegaw=out.nonlinear.signals.values(:,3);
        anglew=out.nonlinear.signals.values(:,4);
    end
    for i=1:length(time)
        
        if ~ishghandle(figCur)
            break
        end
        
        %hold off
        figure(figCur)
        if(get(get(gcf,'currentaxes'),'xlim')~=[-2,2] | get(get(gcf,'currentaxes'),'ylim')~=[-2,2])
            axis([-2,2,-2,2])
        end
        set(plt, 'XData',[0,+sin(fi(i))]);
        set(plt, 'YData',[0,cos(fi(i))]);
        plotCircle(plt,fi(i),anglew(i));
	    pause(0.01);
        
    end
end


function []=plotCircle(plt,fi,anglew)

%plot([sin(fi),sin(fi)+0.5*sin(anglew)],[-cos(fi),-cos(fi)+0.5*cos(anglew)]);
set(plt,'XData',[plt.XData,sin(fi),sin(fi)+0.5*sin(anglew)])
set(plt,'YData',[plt.YData,cos(fi),cos(fi)+0.5*cos(anglew)])
for i=1:36
    %plot([sin(fi)+0.5*sin(anglew+2*pi/36*(i-1)),sin(fi)+0.5*sin(anglew+2*pi/36*i)],[-cos(fi)+0.5*cos(anglew+2*pi/36*(i-1)),-cos(fi)+0.5*cos(anglew+2*pi/36*i)],'r-');
    set(plt,'XData',[plt.XData,sin(fi)+0.5*sin(anglew+2*pi/36*(i-1)),sin(fi)+0.5*sin(anglew+2*pi/36*i)])
    set(plt,'YData',[plt.YData,cos(fi)+0.5*cos(anglew+2*pi/36*(i-1)),cos(fi)+0.5*cos(anglew+2*pi/36*i)])
end

end



function[Nbar]=rscale(a,b,c,d,k)
         % Given the single-input linear system:
         %       .
         %       x = Ax + Bu
         %       y = Cx + Du
         % and the feedback matrix K,
         %
         % the function rscale(sys,K) or rscale(A,B,C,D,K)
         % finds the scale factor N which will
         % eliminate the steady-state error to a step reference
         % for a continuous-time, single-input system
         % with full-state feedback using the schematic below:
         %
         %                         /---------\
         %      R         +     u  | .       |
         %      ---> N --->() ---->| X=Ax+Bu |--> y=Cx ---> y
         %                -|       \---------/
         %                 |             |
         %                 |<---- K <----|
         %
         %8/21/96 Yanjie Sun of the University of Michigan
         %        under the supervision of Prof. D. Tilbury
         %6/12/98 John Yook, Dawn Tilbury revised
         error(nargchk(2,5,nargin));
         % --- Determine which syntax is being used ---
         nargin1 = nargin;
         if (nargin1==2),	% System form
         		[A,B,C,D] = ssdata(a);
         		K=b;
         elseif (nargin1==5), % A,B,C,D matrices
         		A=a; B=b; C=c; D=d; K=k;
         else error('Input must be of the form (sys,K) or (A,B,C,D,K)')
         end;
         % compute Nbar
         s = size(A,1);
         Z = [zeros([1,s]) 1];
         N = inv([A,B;C,D])*Z';
         Nx = N(1:s);
         Nu = N(1+s);
         Nbar=Nu + K*Nx;
end




