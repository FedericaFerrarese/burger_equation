# burger_equation
Finite volume methods and solution to the Burger equation with source term in the one and two dimensional case with different initial data. 

The main.m file solves the Burger equation with source term allowing you to choose among the finite volume methods, defined by flux = num_flux(scheme, lambda, u, v), and the
initial data, defined by z = initialData(x,y,m,dim). The solution at time t is defined by u = sol(u,dt,dx,dy,sc,sp,x,y,dim).

The available numerical methods are upwind, Lax-Wendroff and Rusanov. 
Dependenly on the choice of the initial data, the solution can be smooth or a shock or a rarefraction wave. 
In the one dimensional case, the problem can be solved either fractional method, splitting the equation into its homogeneous part and into an ODE describing the time evolution of the source term. If you decide to use an unsplitted method, the problem will be directly discretized. 
In the two dimensional case, the fractional method is applied not only to split the equation but also to split the dimensions, solving two one-dimensional sub-problems.  
