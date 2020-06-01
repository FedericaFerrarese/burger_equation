%% Burger equation with source 
%
% The script solves the hyperbolic equation 
%                  u_t + f(u)_x = S(x,t)  
% where a = f'(u) is a(u)=u^2, x is a vector, 
% S(x,t)=sum(exp(-norm(x)^2)*dirac(t-tk)) in dimension 1 or 2.
% 
% Methods are defined by flux = num_flux(scheme, lambda, u, v).
% Initial data is defined by z = initialData(x,y,m,dim).
% The solution at time t is defined by u = sol(u,dt,dx,dy,sc,sp,x,y,dim).

clc
clear all
close all

Mx = 60;   % number of x-space steps
My = 80;   % number of y-space steps 
L = 2;     % extrema of the interval [-L,L] (x-axis, y-axis)  

% Select: 
% the dimension of the problem,
% the method to be used,
% if the method is unsplitted or splitted,
% the initial data. 


% One dimension/two dimensions
disp('[1] One dimension');
disp('[2] Two dimensions');

% selection of the dimension of the problem 
dim = input('choose the dimension of the problem (1,2): '); 

% Methods
disp('[1] Upwind');
disp('[2] Lax-Wendroff');
disp('[3] Rusanov');

% selection of the method
sc = input('choose the numerical methods (1-3): '); 

% Unsplitted or splitted 
disp('[1] Unsplitted');
disp('[2] Splitted');

% selection between unsplitted and splitted
sp = input('choose between unsplitted and splitted (1,2): ');

% Set the time and the CFL condition(to ensure stability):
% choose nu such that the CFL condition is satisfied
tf = 3; % final time 
nu = 0.80;   % CFL condition: dt=nu*dx 

% Mesh grid
x = linspace(-L,L,Mx);   
y = linspace(-L,L,My);

% Initial data if the dimension of the problem is 1 
if dim == 1
    disp('[1]  1* (x < 0)+ 0*(x >=0): shock'); % decreasing initial data (shock)
    disp('[2]  -1* (x < 0)+ 1*(x >=0): rarefraction wave'); % increasing initial data (rarefraction)  
    disp('[3] sin(\pi x/L): smooth'); % smooth initial data 
    
    % selection of the initial data 
    m = input('initial data 1-3: ');

    u = initialData(x,y,m,dim);     % intial data 
    a = 2.*u;                       % initial velocity 

    % limitation of the axis
    um  = min(u);                 
    uM  = max(u);
    
    amax = max(abs(a));            % for CFL 
    dx   = (x(end)-x(1)) /(Mx-1);  % x-space step
    dt   = nu*dx/amax;             % time step 
    ntot = tf/dt;                  % total number of timesteps
    
% Initial data if the dimension of the problem is 2     
elseif dim ==2
    disp('[1]  1* (Y < 0)+ 0*(Y >=0): shock'); % decreasing initial data (shock)
    disp('[2]  -1* (Y < 0)+ 1*(Y >=0): rarefraction wave;'); % increasing initial data (rarefraction)  
    disp('[3] sin(\pi (X+Y)/L): smooth;'); % smooth initial data 
    
    % selection of the initial data 
    m = input('initial data 1-3: ');

    u = initialData(x,y,m,dim);     % intial data 
    a = 2.*u;                       % initial velocity 

    % limitation of the axis
    um   = min(u);                
    uM   = max(u);
    um1 = min(um);
    uM1 = max(uM);
    
    amax = max(abs(a));                         % for CFL 
    dx   = (x(end)-x(1)) /(Mx-1);               % x-space step
    dy   = (y(end)-y(1)) /(My-1);               % y-space step 
    dt   = nu*dx*dy./(max(amax).*(dx+dy));      % time step 
    ntot = tf/dt;                               % total number of timesteps
    [X,Y]= meshgrid(y,x);                       % mesh grid
end

tc=0; % initial time


u0 = u; % initial data 

for k = 1:ntot  % time for loop 
   tc = tc+dt;  % time counter
   
   % one dimensional case 
   if dim==1
       dy=0;
       y=0; 
       u  = sol(u,dt,dx,dy,sc,sp,x,y,dim);  % solution at time t=k
   
%      % control CFL condition
%      h=max(2*u);
%      h*dt/dx;
   	
       % Plot in 1 dimension 
       figure(1)
       plot(x,u,':bo');

       s1=sprintf('Solution of the  method %d  at time t=%f',sc,tc);
       title(s1);
       axis([-L L um-.5 uM+.5]);
       pause(.2);
   % two dimensional case     
   elseif dim == 2
       
        u  = sol(u,dt,dx,dy,sc,sp,x,y,dim);  % solution at time t=k  
        
        % control CFL condition
%         amax = max(abs(2*u));
%         h=dt*max(amax)*(dx+dy)/(dx*dy);
        
        % Plot in 2 dimensions  
        figure(1)
        % mesh(X,Y,u)
        %contour(X,Y,u) %level sets
        surf(X,Y,u)
        s1=sprintf('Solution of the method %d at time t=%f',sc,tc);
        xlim([-2 2])
        ylim([-2 2])
        zlim([um1-.5 uM1+.5])
        title(s1);
        pause(.2);
   end
end   
