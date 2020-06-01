% Numerical solution of the Riemann problem at time t, position (x,y)
%
% INPUT:
% u = solution at time t-1
% dt = time step
% dx = x-space step
% dy = y-space step 
% sc = selection of the method:
% 1 = upwind; 2 = Lax-Wendroff; 3 = Rusanov
% sp = selection of splitted or unsplitted method:
% 1 = unsplitted; 2 = splitted 
% x = x-axis discretization
% y = y-axis discretization 
% dim = dimension of the problem 
% 
% OUTPUT: 
% u = solution at time t 

function u=sol(u,dt,dx,dy,sc,sp,x,y,dim)
% one dimensional case  
if dim==1
    lambda= dt/dx; % ratio dt/dx 
    Mx=length(x); % number of x-discretization points 

    % source term 
    source=zeros(1,Mx);
    for j=2:Mx-1
          source(j)=exp(-(x(j).^2));
    end
    r=randi([0 1]);
    s=source*r;
    
    % choose between unsplitted or splitted method 
    if sp == 1 % unsplitted
    % conservative form: u_j^n+1=u_j^n-dt/dx (F_j-1/2^n-F_j-1/2^n) 
        u(2:Mx-1) = u(2:Mx-1)-lambda*(num_flux(sc,lambda,u(2:Mx-1),u(3:Mx))...
            -num_flux(sc,lambda,u(1:Mx-2),u(2:Mx-1)))+dt*s(2:Mx-1);
        
    elseif sp == 2 % splitted 
        % solve without source term 
        u1(2:Mx-1) = u(2:Mx-1)-lambda*(num_flux(sc,lambda,u(2:Mx-1),u(3:Mx))...
            -num_flux(sc,lambda,u(1:Mx-2),u(2:Mx-1))); 
        % introduce the source term 
         u(2:Mx-1) = u1(2:Mx-1) +dt*s(2:Mx-1); % forward-Euler
    end
    
    % Set boundary conditions
    
    % Periodic 
%     u(1)  = u(Mx-1);
%     u(Mx) = u(2);
    % Neumann
    u(1) = u(2);
    u(Mx)= u(Mx-1); 
    
elseif dim==2
    u1=u; % initialization of u1
    u2=u; % initialization of u2 
    lambda=dt/dx; % ratio dt/dx 
    lambda1=dt/dy; % ratio dt/dy
    Mx=length(x); % number of x-discretization points 
    My=length(y); % number of y-discretization points 

    % Source term 
    source=zeros(Mx,My);
    for j=1:My
        for i=1:Mx
             source(i,j)=exp(-(x(i).^2+y(j).^2));
        end
    end
    r=randi([0 1]);
    s=source*r;
    
    % choose between unsplitted or splitted method 
    if sp == 1 % unsplitted 
        for j=2:My-1
            u1(2:Mx-1,j)=u(2:Mx-1,j)-lambda*(num_flux(sc,lambda,u(2:Mx-1,j),u(3:Mx,j))...
                    -num_flux(sc,lambda,u(1:Mx-2,j),u(2:Mx-1,j)))+dt*s(2:Mx-1,j);
        end
        for i=2:Mx-1
            u(i,2:My-1)=u1(i,2:My-1)-lambda1*(num_flux(sc,lambda1,u1(i,2:My-1),u1(i,2+1:My))...
                    -num_flux(sc,lambda1,u1(i,2-1:My-2),u1(i,2:My-1)))+dt*s(i,2:My-1);

        end
    elseif sp==2 % splitted 
         for j=2:My-1 % solution with respect to x-axis
            u1(2:Mx-1,j)=u(2:Mx-1,j)-lambda*(num_flux(sc,lambda,u(2:Mx-1,j),u(3:Mx,j))...
                    -num_flux(sc,lambda,u(1:Mx-2,j),u(2:Mx-1,j)));
        end
        for i=2:Mx-1 % solution with respect to y-axis
            u2(i,2:My-1)=u1(i,2:My-1)-lambda1*(num_flux(sc,lambda1,u1(i,2:My-1),u1(i,2+1:My))...
                    -num_flux(sc,lambda1,u1(i,2-1:My-2),u1(i,2:My-1)));

        end
        % Introduce the source term 
        for j=2:My-1 
             u(2:Mx-1,j) = u2(2:Mx-1,j) +dt*s(2:Mx-1,j);
        end
         for i=2:Mx-1
             u(i,2:My-1) = u2(i,2:My-1) +dt*s(i,2:My-1);
        end
       
    end
    
    % Set boundary conditions
    
    % Periodic 
    u(1,:) = u(Mx-1,:);
    u(:,1) = u(:,My-1); 
    u(Mx,:) = u(2,:);
    u(:,My)=u(:,2);  
    
%    Neumann
%       u(1,:)=u(2,:);
%       u(:,1)=u(:,2);
%       u(Mx,:)=u(Mx-1,:);
%       u(:,My)=u(:,My-1); 
end
end
