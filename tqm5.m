%----------------------------------------------------------------
% tqm5: Calculate the time evolution of a superposition state 
% for a particle trapped in a double-well potential U(x).
%----------------------------------------------------------------
% Parameters
L = 6;
N = 1000;
x = linspace(-L, L, N).';
dx = x(2) - x(1);
% Potential
w = 1*(L/50); % width
a = 1*3*w;  % spacing
U = 0.8*(-100)*( heaviside(x+w-a) - heaviside(x-w-a) ...
         + heaviside(x+w+a) - heaviside(x-w+a)); % potential
% Finite-difference representation of Laplacian and Hamiltonian, 
% where hbar = m = 1.
hbar = 1;
e    = ones(N,1);
Lap  = spdiags([e -2*e e],[-1 0 1],N,N) / dx^2;
H    = -(1/2)*Lap + spdiags(U,0,N,N);
% Find and sort lowest nmodes eigenvectors and eigenvalues of H.
nmodes  = 2; 
[V,E]   = eigs(H,nmodes,'sa'); 
[E,ind] = sort(diag(E)); % Convert E to vector and sort low to high. 
V       = V(:,ind);      % Rearrange coresponding eigenvectors.
% Rescale eigenvectors so that they are always 
% positive at the center of the right well. 
for c = 1:nmodes 
    V(:,c) = V(:,c) / sign(V((3*N/4),c)); 
end
% Compute and display normalized prob. density rho(x,T). 
% Parameters for solving the problem in the interval 0 < T < TF. 
TF = 2*pi*hbar/(E(2)-E(1)); % Length of time interval.
NT = 100;                   % No. of time points.
T  = linspace(0,TF,NT);     % Time vector.
% Compute probability normalization constant (at T=0).
psi_o = 0.5*V(:,1) + 0.5*V(:,2);        % Wavefunction at T=0.
sq_norm = psi_o' * psi_o * dx;          % Square norm = |<ff|ff>|^2.
Usc = 4*U*max(abs(V(:))) / max(abs(U)); % Rescale U for plotting.
% Compute and display rho(x,T) for each time T.
vw = VideoWriter('doubleWell.avi');
open(vw);
for t = 1:NT % time index parameter for stepping through loop.
    % Compute wavefunction psi(x,T) and rho(x,T) at T=T(t).
    psi = 0.5*V(:,1)*exp(-1i*E(1)*T(t)) ...
        + 0.5*V(:,2)*exp(-1i*E(2)*T(t));
    rho = conj(psi) .* psi / sq_norm; % Normalized probability density.
    
    % Plot rho(x,T) and rescaled potential energy Usc.
    plot(x,rho,'-', x, Usc,'-');
    axis([-L/8 L/8 -1 6]);
    xlabel('x (m)');
    ylabel('probability density (1/m)');
    title(['T = ' num2str(T(t), '%05.2f') ' s']);
    
    frame = getframe(gcf);
    writeVideo(vw, frame);
end

close(vw);
