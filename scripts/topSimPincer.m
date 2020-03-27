% topSimPincer.m
%
% Top-level pincer simulation script

%----- Setup
clear;clc;
scenarioVec = {'null','moderate_multipath','severe_multipath', ...
               'matched_power_spoofing','overpowered_spoofing',...
               'weak_jamming','moderate_jamming','severe_jamming'};
scenario = scenarioVec{1};
s.PA = -156;
s.N0 = -204;
s.Ms = 7;
s.taud = 0.15;
s.Ta = 0.1;
s.WFE = 2e6;
s.WFEbeta = 20e6;
s.Pbeta = -130.8815; 
s.Tc = 1e-6;
s.sigmaP = 0.2;

%----- Simulate parameters for scenario
p.N = 1;
p.Tc = s.Tc;
switch scenario
  case 'null'
    p.i = 0;
    op = simulatePincerParameters(p);
    s.Pflag = 0;
  case 'moderate_multipath'
    p.i = 1;
    p.elDeg = 40;
    op = simulatePincerParameters(p);
    s.Pflag = 0;
  case 'severe_multipath'
    p.i = 1;
    p.elDeg = 10;
    op = simulatePincerParameters(p);
    s.Pflag = 0;
  case 'matched_power_spoofing'
    p.i = 2;
    p.elDeg = 5;
    p.muEtadB = 0;
    p.sigmaEtadB = 0.1;
    op = simulatePincerParameters(p);
    s.Pflag = 2;
  case 'overpowered_spoofing'
    p.i = 2;
    p.elDeg = 5;
    p.muEtadB = 10;
    p.sigmaEtadB = 1;
    op = simulatePincerParameters(p);
    s.Pflag = 2;
  case 'weak_jamming'
    p.i = 3;
    p.riceSdB = 3;
    p.riceSigmadB = 1;
    op = simulatePincerParameters(p);
    s.Pflag = 2;
  case 'moderate_jamming'
    p.i = 3;
    p.riceSdB = 10;
    p.riceSigmadB = 1;
    op = simulatePincerParameters(p);
    s.Pflag = 2;
  case 'severe_jamming'
    p.i = 3;
    p.riceSdB = 16;
    p.riceSigmadB = 1;
    op = simulatePincerParameters(p);
    s.Pflag = 2;
end
s.etadB = op.etadB;
s.Delta_tau = op.Delta_tau;
s.Delta_theta = op.Delta_theta;

%----- Simulate observables for scenario
Nsim = 10000;
dVec = zeros(Nsim,1);
PVec = zeros(Nsim,1);
PfullVec = zeros(Nsim,1);
for ii=1:Nsim
  [o] = simulatePincerObservables(s);
  dVec(ii) = o.d; PVec(ii) = o.P; PfullVec(ii) = o.Pfull;
end
% The final symmetric difference D is the modulus of d
DVec = abs(dVec);

%----- Plot results
% Choose P0 so that under scenario "clean" mean(PVec) = P0
P0 = -140.0906; 
fprintf('Scenario: %s\n',scenario);
fprintf('etadB = %.1f dB\n', s.etadB);
fprintf('Delta_tau = %.2f chips\n', s.Delta_tau);
fprintf('Delta_theta = %.2f rad\n', s.Delta_theta);
figure(1);clf;
plot(DVec,PVec - P0, '.');
xlim([0,30]);
ylim([-3,20]);
grid on;
xlabel('D (Normalized IQ units)')
ylabel('P - P_0 (dB)');
title('Scatter plot of excess power and distortion');

figure(2);clf;
subplot(211);
hist(DVec,25);
xlabel('D (Normalized IQ units)');
ylabel('N');
xlim([0 30]);
title('Histograms of distortion and excess power');
subplot(212);
hist(PVec - P0,25);
xlabel('P - P_0 (dB)');
ylabel('N');
xlim([-3 30]);

%----- Check theoretical distribution for D 
% Under the clean or jamming scenarios, the theoretical distribution of D is
% Rayleigh with scale parameter sigma_d^2
if(p.i == 0 || p.i == 3)
  sigma_d = sqrt(8*s.taud*o.beta^2*(o.sigmaN/o.sigmaN0)^2)
  sigma_d_from_samples = std(dVec)
  xVec = [0.001:0.001:7]';
  % Note that the Rayleigh distribution given in Eq. (13) of the paper as posted
  % on IEEE Explore is incorrect.  sigma_d^2 in the equation is twice as large
  % as it should be.  The equation below gives the correct expression.
  pD = (2*xVec/(sigma_d^2)).*exp(-xVec.^2/(sigma_d^2));
  figure(3);clf;
  plot(xVec,pD);
  xlabel('x');
  ylabel('p_D');
  title('Theoretical distribution of D under clean or jamming scenarios');
end