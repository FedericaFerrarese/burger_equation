% Numerical flux functions
%
% INPUT: 
% scheme = selection of the method:
% 1 = upwind; 2 = Lax-Wendroff; 3 = Rusanov
% lambda = dt/dx
% u = u_j
% v = u_j+1
%
% OUTPUT: 
% flux = numerical flux choosen 

function flux = num_flux(scheme, lambda, u, v)

if scheme == 1 % Upwind
    al = u+v;                   % velocity left
    sl = (1+sign(al))/2;        % max(a,0) sl=0 if a<0
    sr = (1-sign(al))/2;        % min(a,0) sr=0 if a>0 
    
    flux = sr.*v.^2+ sl.*u.^2;

    
elseif scheme == 2 % Lax-Wendroff
    flux = 1/2*(u.^2+v.^2)-lambda*(u+v).*(v.^2-u.^2);
    
elseif scheme == 3 % Rusanov
    amax = max(abs(2.*u),abs(2.*v)); 
    flux = 0.5*(u.^2+v.^2)-0.5*amax.*(v-u); 
end   
end
