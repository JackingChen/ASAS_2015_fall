function wavepacket2d(N,dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% wavepacket2d.m
%
% This script simulates the time evolution of a quantum mechanical
% wavepacket on the unit square with Dirichlet boundary conditions.
% Spatial discretization is with a Chebyshev-tau spectral method
% and time discretization is with a second order Crank-Nicholson
% method.
%
% Parameters: 
%
% N - polynomial truncation
% dt - time step size
%
% Written by Greg von Winckel - 09/06/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N1=N+1;
j=0:N;


% Construct derivative matrix
D=zeros(N1,N1);
for k=1:N
    D(k,k+1:2:N1) = 2*(k:2:N);
end
D(1,:)=D(1,:)/2;D2=D*D;
D2=D*D;

% Chebyshev Vandermonde matrix
M=cos(pi*j'*j/N);

% Chebyshev-Gauss-Lobatto nodes
x=cos(pi*j/N);

I=eye(N1);
I2=I; I2(N:N1,:)=0;
II=kron(I2,I2);

% 2D Free particle Hamiltonian (Laplacian)
H=-kron(D2,I2)-kron(I2,D2);
xn=chebnodes(N);

% Chebyshev grid on the square
[x,y]=meshgrid(xn,xn);

% Initial value
u0=exp(-10*((x).^2+(y).^2));

% convert to spectral
U0=fcgltran2d(u0,1);

%%%%% Crank-Nicholson scheme %%%%%
Q=i*H*dt;

% Left and right-hand sides of time stepping operator
LHS=II+Q/2;     RHS=II-Q/2;

% Boundary operator
dir=zeros(N1);
dir(N:N1,:)=M(1:N:N1,:);
BC=kron(dir,I)+kron(I,dir);

% Construct approximate propagator
prop=(LHS+BC)\RHS;

% Inrterpolate on finer mesh using Lagrange polynomials
xf=linspace(-1,1,64)';
[xx,yy]=meshgrid(xf,xf);

% Convert to column
U=U0(:);

for n=1:500

    temp=prop*U;
    U=temp;

    % Convert to point space
    u=fcgltran2d(reshape(U,N1,N1),-1);
    
    % Interpolate data
    u2=barylag2d(abs(u).^2,xn,xn,xf,xf);
    
    % Plot
    mesh(xx,yy,u2);
     axis([-1 1 -1 1 0 1/2]);
     pause(0.001);
end