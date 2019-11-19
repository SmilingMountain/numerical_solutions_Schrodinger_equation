% ----------------------------------------------------
% Implements the shooting method to solve the TISE
% for a particle moving in a 1-D box with V(x) = 0 
% inside domain [0,L] and V(x) infinite elsewhere.
% Atomic units are used.
% Length unit = 1 Angstrom 10^(-10) m
% Mass unit = 1 me (9.1094 x 10^(-31) kg)
% Energy unit = 1 eV (1.6022 x 10^(-19) J)
% In these units the TISE takes the form
% d2u/dx^2 = A*m*(V(x) - E)*u
% Boundary conditions at x = 0 are u = 0, du/dx = 1
% The TISE is written as two 1st order ODE's
% It is integrated using the Matlab solver ode45
% Anonymous function 'rhs' defines the 1sr order DE's
% The unnormalized wavefunction is plotted.
% The trial eigenvalue is set in script
%-----------------------------------------------------
% initial conditions and x range
y0 = [0 1];
xmax = 10;
dx = 0.01;
x = 0:dx:xmax;
% equation parameters
A = 0.26246;
m = 1;
V = 0;
% trial energy eigenvalue for n=1
E_value = 0.37604;
% 1st order representation of the TISE
rhs = @(x,y) [y(2); A*m*(V - E_value).*y(1)];
% call DE solver
[x,y] = ode45(rhs,x,y0);
% normalize solution for u
s = sum(y(:,1).^2);
total = s*dx;
u = y(:,1)/total;
% plot wave function
%plot(x,y(:,1) + E_value);
%hold on;
% plot prob. density
p_d=abs(u).^2;
plot(x,p_d);
hold on;
%plot energy level n=1
E_value1=E_value;
x1=[0 xmax];e1=[E_value1 E_value1];
plot(x1,e1,':');
hold on;
% trial energy eigenvalue for n=2
E_value = 1.504176;
% 1st order representation of the TISE
rhs = @(x,y) [y(2); A*m*(V - E_value).*y(1)];
% call DE solver
[x,y] = ode45(rhs,x,y0);
% normalize solution for u
s = sum(y(:,1).^2);
total = s*dx;
u = y(:,1)/total;
% plot wave function
%plot(x,y(:,1) + E_value);
%hold on;
% plot prob. density
p_d=abs(u).^2;
plot(x,p_d);
hold on;
%plot energy level n=2
E_value2=E_value;
x1=[0 xmax];e1=[E_value2 E_value2];
plot(x1,e1,':');
hold on;
% trial energy eigenvalue for n=3
E_value = 3.38439;
% 1st order representation of the TISE
rhs = @(x,y) [y(2); A*m*(V - E_value).*y(1)];
% call DE solver
[x,y] = ode45(rhs,x,y0);
% normalize solution for u
s = sum(y(:,1).^2);
total = s*dx;
u = y(:,1)/total;
% plot wave function
%plot(x,y(:,1) + E_value);
%hold on;
% plot prob. density
p_d=abs(u).^2;
plot(x,p_d);
hold on;
%plot energy level n=3
E_value3=E_value;
x3=[0 xmax];e3=[E_value3 E_value3];
plot(x3,e3,':');
hold on;
% trial energy eigenvalue for n=4
E_value = 6.016705;
% 1st order representation of the TISE
rhs = @(x,y) [y(2); A*m*(V - E_value).*y(1)];
% call DE solver
[x,y] = ode45(rhs,x,y0);
% normalize solution for u
s = sum(y(:,1).^2);
total = s*dx;
u = y(:,1)/total;
% plot wave function
%plot(x,y(:,1) + E_value);
%hold on;
% plot prob. density
p_d=abs(u).^2;
plot(x,p_d);
hold on;
%plot energy level n=4
E_value4=E_value;
x4=[0 xmax];e4=[E_value4 E_value4];
plot(x4,e4,':');
hold on;
%plot labelling
xlabel('Position in Angstroms');
ylabel('Probability distribution');
%ylabel('psi n');
ttl=['1-D box with shoot1  ']%,...
    %'E1 = ',num2str(E_value1)];
title(ttl);
%axis([0 10 0 0.1])
legend('E1 = ',num2str(E_value1),'E2 = ',num2str(E_value2),'E3 = ',num2str(E_value3),'E4 = ',num2str(E_value4),'location','eastoutside')