function main
%System Set-up %
%Define Variables
%D = 3*10^-8;    % Axial Dispersion coefficient
%v = 1*10^-3;    % Superficial velocity
%epsilon = 0.4;  % Voidage fraction
%k = 3*10^-5;    % Mass Transfer Coefficient
%Kf=2.5*10^-5; %freundlish parameter
%nf= 1.45; %freundlish constant
cFeed = 10;     % Feed concentration
L = 1;         % Column length
t0 = 0;        % Initial Time
%tf = 1000;      % Final time
tf = 2000;      % Final time
dt = 0.5;      % Time step
% z=[0:0.01:L]; %Mesh generation
z = [0:0.0005:L]; %Mesh generation
t = [t0:dt:tf];% Time vector
n = numel(z);  % Size of mesh grid
%Initial Conditions / Vector Creation
c0 = zeros(n,1);     
c0(1) = cFeed;
q0 = zeros(n,1);    % t = 0, q = 0 for all z, this makes sense to me
y0 = [c0 ; q0];     % Appends conditions together
%ODE15S Solver
[T, Y] = ode15s(@(t,y) MyFun(t,y,z,n),t,y0);
%plot(T,Y);
plot(T,Y(:,n)/cFeed)
end 
function DyDt=MyFun(~, y, z, n)
% Defining Constants
D = 3*10^-8;    % Axial Dispersion coefficient
v = 1*10^-3;    % Superficial velocity
epsilon = 0.4;  % Voidage fraction
k = 3*10^-5;    % Mass Transfer Coefficient
Kf=2.5*10^-5; %freundlish parameter
nf= 1.45; %freundlish constant
% Variables being allocated zero vectors
c = zeros(n,1);
q = zeros(n,1);
DcDt = zeros(n,1);
DqDt = zeros(n,1);
DyDt = zeros(2*n,1);
zhalf = zeros(n-1,1);
DcDz = zeros(n,1);
D2cDz2 = zeros(n,1);
c = y(1:n);
q = y(n+1:2*n);
% Interior mesh points
zhalf(1:n-1)=(z(1:n-1)+z(2:n))/2;
for i=2:n-1
  %DcDz(i) = ((z(i)-z(i-1))/(z(i+1)-z(i))*(c(i+1)-c(i))+(z(i+1)-z(i))/(z(i)-z(i-1))*(c(i)-c(i-1)))/(z(i+1)-z(i-1));
  DcDz(i) = (c(i)-c(i-1))/(z(i)-z(i-1));
  D2cDz2(i) = (zhalf(i)*(c(i+1)-c(i))/(z(i+1)-z(i))-zhalf(i-1)*(c(i)-c(i-1))/(z(i)-z(i-1)))/(zhalf(i)-zhalf(i-1));
end
% Calculate DcDz and D2cDz2 at z=L for boundary condition dc/dz = 0
DcDz(n) = 0; 
D2cDz2(n) = -1.0/(z(n)-zhalf(n-1))*(c(n)-c(n-1))/(z(n)-z(n-1));
% Set time derivatives at z=0
% DcDt = 0 since c=cFeed remains constant over time and is already set as initial condition
% Standard setting for q 
DqDt(1) = k*(kf*(c(1)^(1/nf))-q(1));
DcDt(1) = 0.0;
% Set time derivatives in remaining points
for i=2:n
    %Equation: dq/dt = K(q*-q) where q*=kf*(c^(1/nf))
    DqDt(i) = k*(kf*(c(i)^(1/nf))-q(i));
      %Equation: dc/dt = D * d2c/dz2 - v*dc/dz - ((1-e)/(e))*dq/dt
      DcDt(i) = D*D2cDz2(i) - v*DcDz(i) - ((1-epsilon)/(epsilon))*DqDt(i);
  end
% Concatenate vector of time derivatives
DyDt = [DcDt;DqDt];
end
