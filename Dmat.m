%----------------------------------------------------------------
% Numerical calculation of first and second derivatives
% using a differential operator as a matrix.
%----------------------------------------------------------------
% Parameters for solving problem in the interval 0 < x < L.
L = pi;                % Interval Length.
N = 30;                 % No. of coordinate points.
x = linspace(0, L, N).'; % Coordinate vector.
dx = x(2) - x(1);        % Coordinate step.
% Two-point finite-difference representation of Derivative.
D = (diag(ones((N-1),1),1) - diag(ones((N-1),1),-1)) / (2*dx);
% Modify D so that it is consistent with f(0) = f(L) = 0.
D(1,   1) = 0; D(  1, 2) = 0; D(2, 1) = 0; % So that f(0) = 0.
D(N, N-1) = 0; D(N-1, N) = 0; D(N, N) = 0; % So that f(L) = 0.
% Three-point finite-difference representation of Laplacian.
Lap = (diag(ones((N-1),1),-1) - 2*diag(ones(N,1),0) ... 
     + diag(ones((N-1),1), 1)) / (dx^2);
% Modify Lap so that it is consistent with f(0) = f(L) = 0.
Lap(1,   1) = 0; Lap(  1, 2) = 0; Lap(2, 1) = 0; % So that f(0) = 0.
Lap(N, N-1) = 0; Lap(N-1, N) = 0; Lap(N, N) = 0; % So that f(L) = 0.
% To verify that D*f corresponds to taking the derivative of f
% and Lap*f corresponds to taking a second derviative of f,
% define f = sin(x).
f = sin(x); 
% Try the following: 
Df   =   D * f;
Lapf = Lap * f; 
plot(x,f, x,Df, x,Lapf); 
legend('f','Df','Lapf','Location','NorthEast');
xlabel('x');
axis([0 pi -1.1 1.1]);
%daspect([1 1 1]);

% To display the matrix D on screen, simply type D and return ...
disp(D);   % Displays the matrix D in the workspace 
