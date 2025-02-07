%% Dynamics of an aerial ropeway through a time-varying meshing finite element approach

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% INITIAL DATA

%Cabin and cables properties  
Mc= input('Cabin mass [Kg]: '); %[kg] Cabin mass
Ks= input('Carrying cable equivalent stiffness [N/m]: '); %[N/m]Carrying cable equivalent stiffness
T2= input('Tension of hauling cable [N]: '); %[N] Tension of the hauling cable 
m2= input('Mass per unit length of hauling cable [Kg/m]: ');% Mass per unit length of hauling cable [Kg/m]
rho=7800; %[kg/m^3] steel density
A_eff= (m2/rho)*10^6;  %[mm^2] effective section
D_eff= sqrt(4*A_eff/pi); %[mm] effective diameter of the rope

% Geometry and kinematics
v = input('Cablecar velocity [m/s]: '); %[m/s] velocity of the cablecar
L=3;  % [m]ELEMENT LENGTH
t_real= L/v;
q_type_force= input('Which type of force do you want on cabin? (1=inpulse, 2=harmonic):');
if q_type_force==2
force_freq = input('Forcing frequency on cabin [Hz]: '); %[Hz] %force frequency acting vertically on the cabin
end
in_position= input('Cabin initial distance from station [m]: '); % x-coordinate of the cabin node initial position
initial_position= ceil(in_position*201/600); % node of initial position
fin_position= input('Cabin final distance from station [m]: '); % x-coordinate of the cabin node final position
final_position = ceil((fin_position+18)*201/600); % node of final position (18 stays for the meters far from the station in order not to have instability of computation)  
if q_type_force==1
passaggio=0;
p_pilone=input('Insert the position along the span where the input force is happening: ');
pos_pilone=ceil(p_pilone*201/600); % node where there is the impulse force
Mc=0; %the cabin is passing on the pylon, no interia effect of the mass
end
Sim_anim= input('Do you want to see the animation (1) or simulation analysis/plots (2): ');
if Sim_anim==1
q_save= input('Do you wan to save the animation? 1=Yes, 2=No: ');
else
q_save=2;
end
if Sim_anim==2
q_traj=input('Do you want to see the vertical trajectory of the cabin? (1=yes, 2=no):');
q_mids=input('Do you want to see the vertical trajectory of the cable midspan? (1=yes, 2=no):');
q_stress= input('Do you want to monitor Von mises stresses? (1=Yes, 2=No): ');
else 
q_stress=2;
end

% Time scale definition
time_scale= v/(in_position*((sqrt(T2/m2)/(2*in_position))));  % time scale separation 
dt=t_real*time_scale; % time discretization
step_position=1; % Variation of cable length between two successive time simulations
time_to_step= step_position*L/v;  % time for the cabin to travel step_position

% Check on the maximum lenght of the hauling cable element
wmax= input('Maximum frequency of the analysisn [Hz]: ');
eta= 1.5; % safety coefficient on the maximum allowable length 
Lmax2 = (3.14/(eta*wmax))*(T2/m2)^0.5; % maximum length of a single element of the hauling cable
disp(['Maximum allowable cable element length [m]: ' num2str(Lmax2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD INITIAL GEOMETRY AND ASSEMBLE THE SYSTEM
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,T,posit,nbeam,pr]=loadstructure;
indice=0;
CONTATORE=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ASSIGN INITIAL CONDITIONS
nnod_i0= initial_position+1; %number of nodes at the initial time istant (+1 due to the adding of the node for the spring between the cabin and the carrying cable)
ndof_i0=nnod_i0*3-6; % number of  degree of freedon at initial time instant
x0 = zeros(ndof_i0,1); % vertical displacement of the haul rope is 0 at first time istant
xdot0= zeros(ndof_i0,1); % vertical velocity of the haul rope is 0 at forst time istant
kk=0;
gamm = 1/2; %gamma coefficient fot Newmark integration function
beta = 1/4; %beta coefficient fot Newmark integration function
trajectory=[]; %initializing the trajectory vector
midspan=[]; %initailising the midspan vector

if Sim_anim==1
% Animation parameters
num_frames = floor((((initial_position*L)-((final_position-1)*L))/v)/dt); % Numero di frame dell'animazione
fig = figure;

% Iniatlizinga array to create the animation
frames(num_frames) = struct('cdata', [], 'colormap', []);

jjj=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIMULATION
%initialising vectors for simulation
indice=0;
CONTATORE=0;
Y=[];
O=[];
Moment=[];
Shear=[];
Cast_socket_moment=[];
Cast_socket_shear=[];

for pos_cab= initial_position:-step_position:final_position %change of geometry every step_position
%assigning the temporary values to all input for assem function (function that created temporary M and K matrix)
xy_i= xy(1:1:pos_cab,:); 
xy_i(pos_cab+1,1)=xy_i(pos_cab,1); %additional node for the spring
xy_i(pos_cab+1,2)=200; %additional node for the spring
nnod_i= pos_cab+1;
sizee_i= sizee;
ndof_i=nnod_i*3-6;
l_i=l(1:1:(pos_cab-1));
gamma_i=gamma(1:1:(pos_cab-1));
m_i=m(1:1:(pos_cab-1));
EA_i=EA(1:1:(pos_cab-1));
EJ_i=EJ(1:1:(pos_cab-1));
T_i=T(1:1:(pos_cab-1));
posit_i= posit(1:1:(pos_cab-1),:);
nbeam_i= pos_cab-1;
pr_i=pr(1:1:(pos_cab-1));

%building idb matrix (NUMBERS ARE ASSIGNED PROGRESSIVELY FIRST TO FREE
%NODES THEN PROGRESSIVELY TO CONSTRAINED NODES)
idb_i=zeros(nnod_i,3); %initializing idb matrix
%free nodes
idb_i(1,3)=idb(1,3); % free rotation of the first node of traction cable
idb_i(2:1:(nnod_i-2),:)=idb(2:1:(nnod_i-2),:); % free nodes of traction cable
idb_i(nnod_i-1,2)=idb_i(nnod_i-2,3)+1; %free y displacement of the last node of traction cable
idb_i(nnod_i-1,3)=idb_i(nnod_i-2,3)+2; %free rotation of the last node of traction cable
%constrained nodes (-6 dof)
idb_i(1,1)=idb_i(nnod_i-2,3)+3; %constrained x movement of the first node
idb_i(1,2)=idb_i(nnod_i-2,3)+4; %constrained y movement of the first node
idb_i(nnod_i-1,1)=idb_i(nnod_i-2,3)+5; %constrained x movement of the cabin
idb_i(nnod_i,1)=idb_i(nnod_i-2,3)+6; %constrained x movement of the cabin spring
idb_i(nnod_i,2)=idb_i(nnod_i-2,3)+7; %constrained y movement of the cabin spring
idb_i(nnod_i,3)=idb_i(nnod_i-2,3)+8; %constrained rot movement of the cabin spring

%building incid matrix
incid_i=zeros(nbeam_i,6);
for j=1:1:(nbeam_i)  
    incid_i(j,4:1:6)=idb_i(j,:);
    incid_i(j,1:1:3)=idb_i(j+1,:);
end

j=0;
%% Draw structure
%dis_stru(posit_i,l_i,gamma_i,xy_i,pr_i,idb_i,ndof_i);

%% Assemble mass  stiffness and damping matricies
[M,K] = assem(incid_i,l_i,m_i,EA_i,EJ_i,T_i,gamma_i,idb_i,v);

%Add concentrated elements
%  spring 
i_ndof_spring_t = idb_i(pos_cab,2);
i_ndof_spring_s = idb_i(pos_cab+1,2);
i_dofK = [i_ndof_spring_s i_ndof_spring_t];
K_int = [Ks -Ks; -Ks  Ks];
K(i_dofK,i_dofK) = K(i_dofK,i_dofK) + K_int;

% mass
i_ndof2= idb_i(pos_cab,2);
M(i_ndof2,i_ndof2) = M(i_ndof2,i_ndof2) + Mc;

% Compute natural frequencies and mode shapes
MFF = M(1:ndof_i,1:ndof_i);
MCF = M(ndof_i+1:end,1:ndof_i);
MFC = M(1:ndof_i,ndof_i+1:end);
MCC = M(ndof_i+1:end,ndof_i+1:end);

KFF = K(1:ndof_i,1:ndof_i);
KCF = K(ndof_i+1:end,1:ndof_i);
KFC = K(1:ndof_i,ndof_i+1:end);
KCC = K(ndof_i+1:end,ndof_i+1:end);

%building damping matrix C (Rayleigh hypothesys)
ab= [0.008 0.008];  %Rayleigh coefficients
C = ab(1)*M + ab(2)*K;
CFF = C(1:ndof_i,1:ndof_i);
CCF = C(ndof_i+1:end,1:ndof_i);
CFC = C(1:ndof_i,ndof_i+1:end);
CCC = C(ndof_i+1:end,ndof_i+1:end);

%% Force addition

if q_type_force==2  %harmonic force
simulation_time = time_to_step;  %[s] the single simulation time is equal to the time neede to travel one step_position
t = 0:dt:simulation_time; %time discretization 
R = zeros(ndof_i,length(t)); %initializing force vector
jj=1; %time index on the force vector R
for ii=0:dt:simulation_time
R(idb_i(pos_cab,2),jj)=300*sin(2*pi*force_freq*dt*kk); %harmonic force
kk=kk+1; %index of time along the force (it doesn't go to 1 at the end of each single simulation to keep the phase consistent between to successive simulations)
jj=jj+1;
end
end

if q_type_force==1  %inpulse force
simulation_time = time_to_step;  %[s] the single simulation time is equal to the time neede to travel one step_position
t = 0:dt:simulation_time; %time discretization 
R = zeros(ndof_i,length(t)); %initializing force vector
if pos_cab<=pos_pilone 
if passaggio==0
R(idb_i(pos_cab,2),2)=500; %IMPULSE FORCE
R(idb_i(pos_cab,2),3)=1000;
R(idb_i(pos_cab,2),4)=2000;
R(idb_i(pos_cab,2),5)=1000;
R(idb_i(pos_cab,2),6)=1000;
R(idb_i(pos_cab,2),7)=800;
R(idb_i(pos_cab,2),8)=600;
R(idb_i(pos_cab,2),9)=400;
R(idb_i(pos_cab,2),10)=200;
passaggio=1;
end
end
end


%% Time integration (Newmark method)
[x,dx,ddx] = NewmarkSolverMultiple(MFF,CFF,KFF,R,x0,xdot0,t,dt,gamm,beta);

%consistency between two successive simulations:
x0=x((3*step_position+1):1:ndof_i,end);   % inital displacemet for simulation n+1 = final displacement of simulation n - displacement of the nodes passed inside the machine room during a single simulation
xdot0=dx((3*step_position+1):1:ndof_i,end); % inital velocities for simulation n+1 = final velocities of simulation n - velocities of the nodes passed inside the machine room during a single simulation
%%
%Weighted mean
x0(1)=(x(4,end)+4*x(1,end))/5;
x0(2)=(x(5,end)+4*x(2,end))/5;
x0(3)=(x(6,end)+4*x(3,end))/5;
xdot0(1)=(dx(4,end)+4*dx(1,end))/5;
xdot0(2)=(dx(5,end)+4*dx(2,end))/5;
xdot0(3)=(dx(6,end)+4*dx(3,end))/5;


x0(4)=(2*x(7,end)+3*x(4,end))/5;
x0(5)=(2*x(8,end)+3*x(5,end))/5;
x0(6)=(2*x(9,end)+3*x(6,end))/5;
xdot0(4)=(2*dx(7,end)+3*dx(4,end))/5;
xdot0(5)=(2*dx(8,end)+3*dx(5,end))/5;
xdot0(6)=(2*dx(9,end)+3*dx(6,end))/5;

x0(7)=(3*x(10,end)+2*x(7,end))/5;
x0(8)=(3*x(11,end)+2*x(8,end))/5;
x0(9)=(3*x(12,end)+2*x(9,end))/5;
xdot0(7)=(3*dx(10,end)+2*dx(7,end))/5;
xdot0(8)=(3*dx(11,end)+2*dx(8,end))/5;
xdot0(9)=(3*dx(12,end)+2*dx(9,end))/5;

x0(10)=(4*x(13,end)+1*x(10,end))/5;
x0(11)=(4*x(14,end)+1*x(11,end))/5;
x0(12)=(4*x(15,end)+1*x(12,end))/5;
xdot0(10)=(4*dx(13,end)+1*dx(10,end))/5;
xdot0(11)=(4*dx(14,end)+1*dx(11,end))/5;
xdot0(12)=(4*dx(15,end)+1*dx(12,end))/5;

%% Calculation of M and T
if q_stress==1
Y=x(idb_i(2:end-1,2),:);   %saving the vertical displacement history of the cable nodes 
O=x(idb_i(2:end-1,3),:);   % saving the rotation displacement history of the cable nodes 
[mm,nn]=size(Y);

for j=1:nn
for i=1:mm-1
    c=-(3/(L*L))*Y(i,j)+(3/(L*L))*Y(i+1,j)-(2/L)*O(i,j)-(1/L)*O(i+1,j);
    d=(2/(L*L*L))*Y(i,j)-(2/(L*L*L))*Y(i+1,j)+(1/(L*L))*O(i,j)+(1/(L*L))*O(i+1,j);
    if abs(2*c+6*d*L)>abs(2*c)         %calculation of the moment time history for each beam
    Moment(i,j)=EJ_i(1)*(2*c+6*d*L);  
    else
    Moment(i,j)=EJ_i(1)*(2*c);
    end
    Shear(i,j)=EJ_i(1)*6*d; %calculation of the shear time history for each beam
end
end


%% Searching beam that has the maximum sollecitations  (it is always the last one, the nearest to the cabin (where there is the cast socket))
[Max_moment, beam_m ]= max(Moment); % maximum sollecitation and relative number of the beam
[Max_shear, beam_s ]= max(Shear); % maximum sollecitation and relative number of the beam
Cast_socket_moment = [Cast_socket_moment,  Moment(end,:)]; %collecting the total time history of the cast socket segment moment
Cast_socket_shear =  [Cast_socket_shear, Shear(end,:)]; %collecting the total time history of the cast socket segment shear


Moment=[];
Shear=[];
if q_stress==1
sigma_a= T(1)/A_eff; %[MPa] stress due to axial tension
 sigma_m= (32*1000*Cast_socket_moment)/(pi*D_eff^3); %[Mpa] stress due to bending in the most stressed point
 max_sigma_m= max(abs(sigma_m)); %[Mpa] maximum stress due to bending 
 sigma_am= sigma_a + sigma_m; %[Mpa] stress due to bending amd axial tension in the most stressed point
 tau= 16*(Cast_socket_shear)/(3*pi*D_eff^2); %[Mpa]  stress due to shear in the most stressed point
 sigma_VM= sqrt((sigma_am.^2)+3*(tau.^2)); %[MPa] Von Mises stress in the most stressed point
 max_sigma_VM= max(sigma_VM); %[MPa] Max Von Mises stress in the most stressed point
end
end

trajectory= [trajectory, x(idb_i(pos_cab,2),(1:1:end-1))]; % building the  vertical trajecotry history of the cabin 
midspan=[midspan, x(idb_i(ceil(pos_cab/2),2),(1:1:end-1))]; % building the  vertical trajecotry history of a point of the rope a the middle between the cain and the station
CONTATORE=CONTATORE+1;

if Sim_anim==1
for i = 1:((L/v)/dt)
    % Creation of the plot for the frame
    plot(xy_i(2:1:(end-1),1),x(idb_i(2:1:(end-1),2),i))
    rectangle ('Position',[xy_i(end-1,1) - 10, x(idb_i(pos_cab,2),i) - 0.1, 20,0.07], 'FaceColor','black');
    xlim([0 600])
    ylim([-1 1])
    xlabel('L(t) [m]')
    ylabel('[m]')
    title('Entire system vibrations')
    % Save the frame in frame array
    frames(jjj) = getframe(fig);
    jjj=jjj+1;
    % Clean plot to get the new file
    clf;
end
end
Computation_percentage= 100-(pos_cab-final_position)/(initial_position-final_position)*100  %percentage of computation
end

disp('Simulation computed successfully!')
t_sim= (CONTATORE-1)*step_position*L/v;  % [s] total time to travel from initial position to final position
time= linspace(0,t_sim,length(trajectory)); % discretization of the total time
space=linspace((final_position-1)*L,(initial_position-1)*L,length(trajectory)); %discretization of the total space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT
if Sim_anim==2

if q_traj==1
figure()
plot(space,flip(trajectory))
xlabel('L(t) [m]', 'FontSize', 14)
ylabel('y cabin [m]', 'FontSize', 14)
xlim([0 600])
ylim([-1 1])
title(['Trajectory of the cabin, with v=',num2str(v),' m/s'])
end

if q_mids==1
figure()
plot(space/2,flip(midspan))
xlabel('L(t) [m]', 'FontSize', 14)
ylabel('y midspan [m]', 'FontSize', 14)
xlim([0 300])
ylim([-1.5 1.5])
title(['Trajectory of midspan, with v=',num2str(v),' m/s'])
end

if q_stress==1 
figure()
time=linspace(0,t_sim,length(sigma_VM));
plot(time,sigma_VM)
xlabel('Simulation time [s]', 'FontSize', 10)
ylabel('Von Mises Stress [MPa]', 'FontSize', 10)
title(['Von Mises Stress,v=',num2str(v),' m/s'])
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q_save==1

    % Create a file video 
videoFile = VideoWriter('animazione_video.avi', 'Uncompressed AVI');
videoFile.FrameRate = num_frames/t_sim ; % Imposta la frequenza di frame desiderata (frame al secondo)
open(videoFile);

% Write each frame in the animation
for i = 1:num_frames
    writeVideo(videoFile, frames(i).cdata);
end

% Close file video 
close(videoFile);

% Close plot window 
close(fig);


end






