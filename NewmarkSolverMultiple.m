function [u,u_1p,u_2p] = NewmarkSolverMultiple(M,C,K,P,u0,u0_1p,tv,dt,gamma,beta)
%==========================================================================
% Newmark's Direct Integration Method
%==========================================================================
% Solves a multiple-DOF system with mass matrix "M", spring stiffness 
% matrix "K" and damping coeffiecient matrix "C", under exciting force 
% matrix "P".
%
% Returns the displacement, velocity and acceleration of the system
%==========================================================================
% SYNTAX
%       [u,u_1p,u_2p] = NewmarkSolver(M,C,K,P,u0,u0_1p,tv,dt,gamma,beta)
%==========================================================================
% INPUTs
%       [M]:        System Mass              [length(M),length(M)]
%       [C]:        System Damping           [length(M),length(M)]
%       [K]:        System Stiffness         [length(M),length(M)]
%       [P]:        Externally Applied Load  [length(M),length(tv)]
%       [u0]:       Initial Position         [length(M),1]
%       [u0_1p]:    Initial Velocity         [length(M),1]
%       [tv]:       Time Vector              [length(tv),1]
%       [dt]:       Time step                Scalar
%       [gamma]:    Parameter                Scalar
%       [beta]:     Parameter                Scalar
%==========================================================================
% OUTPUTs
%       [u]:        Displacemente Response   [length(M),length(tv)]
%       [u_1p]:     Velocity                 [length(M),length(tv)]
%       [u_2p]:     Acceleration             [length(M),length(tv)]
%==========================================================================
%       By: Xingzhuang Zhao (timezhaox@gmail.com) 2020.12
%==========================================================================
% NOTE:
% Please use gamma = 1/2 and beta = 1/4 as the key parameters or refer to 
% Chopra, A. K. (2020), Dynamics of structures (Fifth Edition), p168 for
% more information about Newmark coefficients and time step-size
%==========================================================================

% Initial setup
u = zeros(length(K),length(tv));
u_1p = zeros(length(K),length(tv));
u_2p = zeros(length(K),length(tv));

% Set the intial conditions
u(:,1) = u0;
u_1p(:,1) = u0_1p;
u_2p(:,1) = M\(P(:,1)-C*u_1p(:,1)-K*u(:,1));

% Computing constants
a1 = M/(beta*dt^2) +  C*gamma/(beta*dt);
a2 = M/(beta*dt) + (gamma/beta-1)*C;
a3 = M*(1/(2*beta)-1) + C*dt*(gamma/(2*beta)-1);
kb = K + a1;

% Step-by-step calculations
for i = 1:length(tv)-1
    pb = P(:,i+1) + a1*u(:,i) + a2*u_1p(:,i) + a3*u_2p(:,i);
    u(:,i+1) = kb\pb;
    u_1p(:,i+1) = gamma/(beta*dt)*(u(:,i+1)-u(:,i)) + (1-gamma/beta)*u_1p(:,i) + dt*(1-gamma/(2*beta))*u_2p(:,i);
    u_2p(:,i+1) = (u(:,i+1)-u(:,i))/(beta*dt^2) - u_1p(:,i)/(beta*dt) - (1/(2*beta)-1)*u_2p(:,i);
end

end
