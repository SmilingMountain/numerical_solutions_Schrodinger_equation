%heaviside example
L = 5;                   % Interval Length.
N = 1000;                % No of points.
x = linspace(-L, L, N).';% Coordinate vector.
dx = x(2) - x(1);        % Coordinate step.
% Finite square well of width 2w and depth given.
w = L/10;
U = -1*(heaviside(x+w)-heaviside(x-w));
figure(1);
plot(x,U);
hold on;
U2 = -1*heaviside(x+w)+1.51;
U3 = -1*heaviside(x-w)+1.49;
plot(x,U2,'r',x,U3,'g');
axis([-L L -2 2]);
hold off;
figure(2);
%Two finite square wells of width 2w and distance 2a apart.
a = 6*w;
U4 = -1*(heaviside(x+w-a) - heaviside(x-w-a) ... 
  + heaviside(x+w+a) - heaviside(x-w+a));
plot(x,U4);
axis([-L L -2 2]);
