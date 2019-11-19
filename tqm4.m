%----------------------------------------------------------------
% tqm4: Finds the first few lowest bound eigenfunctions psi(x)
% and eigenvalues E of 1D TISE,
% -1/2*hbar^2/m*(d2/dx2)psi(x) + U(x)*psi(x) = E*psi(x),
% for a given analytic potential U(x).
% Units are chosen so that hbar = m = 1.
%----------------------------------------------------------------
% Parameters for solving problem in the interval -L < x < L.
L = 5;                      % Interval Length.
N = 1000;                   % No of points.
x = linspace(-L, L, N).';   % Coordinate vector.
dx = x(2) - x(1);           % Coordinate step.
% Potential: 4 options are given below
%U = 1/2*100*x.^(2); % quadratic harmonic oscillator potential
a=0;
b=1.5;
U = (10/(b^4))*((x-(a/2)).^2-b^2).^2-100; % quartic potential
% Finite square well of width 2w and depth given.
%w = L/5;
%U = -10*(heaviside(x+w)-heaviside(x-w));
% Two finite square wells of width 2w and distance 2a apart.
%w = L/20;
%a = 0.6*1.5*w;
%U = -100*(heaviside(x+w-a) - heaviside(x-w-a) ... 
% + heaviside(x+w+a) - heaviside(x-w+a));
% Three-point finite-difference representation of Laplacian
% using sparse matrices that save memory by only
% storing non-zero matrix elements.
e = ones(N,1);
Lap = spdiags([e -2*e e],[-1 0 1],N,N) / dx^2;
% Total Hamiltonian.
hbar = 1;
m = 1;
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
% Find lowest nmodes eigenvectors and eigenvalues of sparse matrix.
nmodes = 3; 
[V,E] = eigs(H,nmodes,'sa'); % find eigenvalues and eigenvectors.
[E,ind] = sort(diag(E));     % convert E to vector and sort low to high.
V = V(:,ind);                % rearrange corresponding eigenvectors.
% Generate plot of lowest energy eigenvectors V(x) and U(x).
Usc = 10 * U * max(abs(V(:))) / max(abs(U));% rescale U for plotting.
plot(x,-V, x,Usc,'--k');     % plot V(x) and rescaled U(x).
xlabel('x');
ylabel('psi(x)');
title('unnormalized wavefunctions of the displayed potential');
% Add legend showing Energy of plotted V(x).
legendLabels = [repmat('E = ',nmodes,1), num2str(E)];
legend(legendLabels) % place lengend string on plot.
xlim([-5 5]);
ylim([-0.2 0.2]);

