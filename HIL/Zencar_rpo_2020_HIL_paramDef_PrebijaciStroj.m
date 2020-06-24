clc
clear all
close all

g=9.81;


%reloading machine
initial_pos=0.25;%must be between 0 and stroke(0.3)
m_sys=5;
stroke=0.3;
b=0.0002;

%gearbox
GR_rot=1/53;
GR_transl=0.05;%rack and pinion wheel radius %1 revolution = 31.42cm (2*pi*GR_transl)

%encoder
rot_enc_ppr=2048;

%HAL sensor
hal_maxi=100;
hal_maxv=5;
hal_zerov=2.5;

%motor
Umax=24;
Ra=0.102;
La=0.016e-3;
cfi=0.0137;
J_motor=3.3e-6;
J_gearbox=9.0000e-07;
J_endgear_and_press=m_sys*GR_transl^2*GR_rot^2;
J=J_motor+J_gearbox+J_endgear_and_press;

GR_overall=GR_rot*GR_transl;

%DC motor state space
A=[0 1 0;0 -b/J cfi/J;0 -cfi/La -Ra/La];
B=[0;0;1/La];
C=eye(3);
D=[0;0;0];

Q=diag([50,5e-4,5e-4]);
R=1;
K=lqr(A,B,Q,R);

Nbar=rscale(A,B,[1 0 0],0,K);







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


