clear all
% Data ====================================================================
eta = 0.78;
Dp = 0.001;
Dz = 9.89019E-05;
Deff = 0.001; %unused
u = 0.036;
C0 = 0.0437; %mg/L
rho = 997.3;
KL = 0.63;  %1/mg
qmax = 154; %mg/g
kf = 0.027419644;
rho_b = 743;


% Axial and radial grid points ============================================
L = 0.04; %m
AxialGridPoints = 201;
InteriorAxialGridPoints = AxialGridPoints - 1;
dz = L/AxialGridPoints;

rp = 1.1/1000; %m
ParticleGridPoints = 21;
InteriorParticleGridPoints = ParticleGridPoints - 1;
dr = rp/ParticleGridPoints;
ParticleCoordinate = linspace(0,rp,ParticleGridPoints);


% Time ====================================================================
t = linspace(0,1440,100);


% Initial and Boundary conditions =========================================
% reactor
C = zeros(AxialGridPoints,1);
C(1) = C0; %z = 0
C(end) = C(end-1); %z = L 

% Particle
q = zeros(AxialGridPoints, ParticleGridPoints);
q(:,1) = q(:,2); %r = 0
coefficient = (2*kf*dr)/(Dp);
q(:,ParticleGridPoints) = (coefficient.*C + 4.*q(:,InteriorParticleGridPoints) -...
    q(:,InteriorParticleGridPoints-1))./(3 + coefficient); %Second order backwards at r = rp 
                                                         % (Might be unnecessary just yet)


% Packing data ============================================================
Initial = [C q];
Data = [eta,Dp,Dz,u,C0,rho,KL,qmax,kf,rho_b];
Spatial = [AxialGridPoints,InteriorAxialGridPoints,ParticleGridPoints,... 
           InteriorParticleGridPoints,ParticleCoordinate,dz,dr,rp];


% ODE solver ==============================================================

[t, y] = ode15s(@PoreDiffusion,t,Initial,[],Data,Spatial);
beep

% Result ==================================================================
CFinal = y(:,1:AxialGridPoints);
qFinal = y(:,AxialGridPoints+1:end); 

% Loop ====================================================================

function dydt = PoreDiffusion(t,y,Data,Spatial)
        
        %Unpacking data---------------------------------------------------- 
        eta = Data(1);
        Dp = Data(2);
        Dz = Data(3);
        u = Data(4);
        C0 = Data(5); 
        rho = Data(6);
        KL = Data(7);
        qmax = Data(8);
        kf = Data(9);
        rho_b = Data(10);


        AxialGridPoints = Spatial(1);
        InteriorAxialGridPoints = Spatial(2);
        ParticleGridPoints = Spatial(3);
        InteriorParticleGridPoints = Spatial(4);
        ParticleCoordinate = Spatial(5:length(Spatial)-3);
        dz = Spatial(length(Spatial)-2);
        dr = Spatial(length(Spatial)-1);
        rp = Spatial(length(Spatial));
        %------------------------------------------------------------------

        % Sorting out input -----------------------------------------------
        C = y(1:AxialGridPoints); % 1 - 201
        q = y(AxialGridPoints+1:end); % 202 - end
        q = reshape(q,AxialGridPoints,ParticleGridPoints);
        M = (3*(1-eta)*kf)/(rp*eta);
        coefficient = (2*kf*dr)/(Dp);

        Cs = zeros(AxialGridPoints,1);
        dqdt = zeros(AxialGridPoints,ParticleGridPoints);
        dCdt = zeros(AxialGridPoints,1);
        % -----------------------------------------------------------------

        for z = 2:InteriorAxialGridPoints % z is the axial coordinate

            % Particle
            for r = 2:InteriorParticleGridPoints
                dqdr = q(z,r)/2/dr - q(z,r)/2/dr;
                d2qdr2 = (q(z,r+1) - 2*q(z,r) + q(z,r-1))/dr/dr;
                dqdt(z,r) = (Dp/ParticleCoordinate(r))*(2*ParticleCoordinate(r)*dqdr ...
                    + ParticleCoordinate(r)^2*d2qdr2);
            end

            % Reactor
            dqdt(:,end) = (coefficient.*C + 4.*q(:,InteriorParticleGridPoints) -...
                            q(:,InteriorParticleGridPoints-1))./(3 + coefficient); %Second order backwards
            Cs(z) = dqdt(z,end)/(KL*qmax - KL*dqdt(z,end)); % q instead of dqdt?
            dCdz = C(z+1)/2/dz - C(z-1)/2/dz;
            d2Cdz2 = (C(z+1) - 2*C(z) + C(z-1))/dz/dz;
            dCdt(z) = -(u/eta)*dCdz + (Dz/eta)*d2Cdz2 - M.*(C(z) - Cs(z));
        end
        
      
        % Reassigning boundary conditions
        dCdt(AxialGridPoints) = dCdt(InteriorAxialGridPoints); %z = L
        dqdt(:,1) = Dp.*(q(:,3) - 2.*q(:,2) + q(:,1))/dr^2; %r = 0
        
        
        dqdt = reshape(dqdt,[],1);
        dydt = [dCdt;dqdt];
end


