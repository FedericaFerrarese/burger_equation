% Intial data for hyperbolic problems
%
% INPUT:
% x = x-space interval
% y = y-space interval 
% m = selection of the initial data
% dim = dimension of the problem
% 
% OUTPUT: 
% z = initial data choosen 

function z = initialData(x,y,m,dim)

% one dimensional case  
if dim==1 
    L=max(x); % extrema of the space interval 
    n=length(x); % number of space steps 
    z=zeros(n); % initialization of the initial data 
    
    if m==1
        z = 1 * (x < 0)+  0* (x >=0);  % shock case  
    elseif m==2
        z = -1* (x < 0)+ 1*(x >=0);  % rarefraction wave case 
    elseif m==3
         z = sin(pi*x/L); % smooth case  
    end
    
% two dimensional case 
elseif dim==2  
    L=max(x); % extrema of the mesh grid 
    [X,Y]=meshgrid(y,x); % mesh grid (spatial discretization) 
    
    if m==1
        z = 1 * (Y < 0)+0 * (Y >=0);  % shock case  
    elseif m==2
        z = -1* (Y < 0)+ 1*(Y >=0);  % rarefraction wave case 
    elseif m==3
         z = sin(pi*(X+Y)/L); % smooth case  
    end
end