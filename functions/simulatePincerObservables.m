function [o] = simulatePincerObservables(s)
% simulatePincerObservables : Simulate the observables used for hypothesis
%                             testing in the pincer defense.
%
%
% INPUTS
%
% s ------- structure variable with the following fields:
% 
%               PA = Power in authentic signal within bandwidth WFEbeta (dBW)
%
%            etadB = Power advantage of interference signal (dB)
%
%        Delta_tau = tau_I - tau_A, the difference between the interference
%                    and the authentic code phases (chips)
%
%      Delta_theta = theta_I - theta_A, the difference between the
%                    interference and the authentic carrier phase (rad)
%
%               N0 = Thermal noise floor (dBW/Hz)
%
%                M = Number of multi-access signals (excluding the
%                    authentic desired signal) with average power PA
%      
%             taud = Symmetric difference measurement offset from tau = 0
%                    (chips).  
%
%               Ta = Accumulation (averaging) interval (seconds)
%
%              WFE = Bandwidth over which received in-band power is measured
%                    (Hz)
%
%          WFEbeta = Bandwidth over which the AGC operates in calculating the
%                    scaling factor beta.  WFEbeta >= WFE may be equal to WFE,
%                    but doesn't have to be.  Typically, WFEbeta is equal to
%                    the full complex bandwidth of the front end (e.g., 20 MHz
%                    for the 25 Msps TEXBAT data).
%
%            Pbeta = Setpoint for AGC: The incoming signal r will be scaled
%                    by beta such that the power in r*beta is Pbeta (dBW)
%
%               Tc = Chip interval (seconds)
%
%           sigmaP = Standard deviation of in-band power measurement.  This
%                    deviation results from measurement error due to
%                    thermal noise, variation in power due to
%                    constellation movement, and variation due to
%                    temperature changes in the receiver LNA (dB). 
%
%            Pflag = Indicates whether the average power in the (M-1)
%                    multi-access signals (excluding the desired signal) is
%                    taken to be equal to the authentic power PA (Pflag = 0),
%                    or equal to the power of the coherently-combined
%                    authentic and interference signals PIA (Pflag = 1), or
%                    equal to the non-coherently-combined authentic and
%                    interference power (Pflag = 2).  Set Pflag = 1 for a
%                    phase-coherent spoofing attack in which all channels are
%                    subject to attack with the same relative phase offset.
%                    Set Pflag = 2 for a spoofing attack or multipath
%                    situation in which the multi-access signals (excluding
%                    the desired signal) are subject to a spoofing (or
%                    multipath) signal having the same interference power as
%                    the interference signal corresponding to the desired
%                    signal, but with (uniformly) random relative phasing.
%                    Set Pflag = 0 for a spoofing attack against a single
%                    signal, or a multipath situation in which it is assumed
%                    that only the desired signal is experiencing significant
%                    multipath.
%
%
% OUTPUTS
%
% o ------- structure variable with the following fields:
%
%        d = Simulated complex symmetric difference scaled by 1/sigmaN0 
%
%        P = Simulated received power measurement in the band WFE (dBW)
%
%    Pfull = Total received power in the band WFEbeta >= WFE (dBW)
%
%+------------------------------------------------------------------------------+
% References: 
%
% Author:  Todd Humphreys
%+==============================================================================+


%----- Calculate beta according to the received power model
eta = 10^(s.etadB/10);
PAL = 10^(s.PA/10); 
N0L = 10^(s.N0/10);
PIAL = (1+eta)*PAL + 2*sqrt(eta)*PAL*cos(s.Delta_theta)*Rcorr(s.Delta_tau);
if(s.Pflag == 0)
  PbarMultiAccess = PAL;
elseif(s.Pflag == 1)
  PbarMultiAccess = PIAL;
elseif(s.Pflag == 2)
  PbarMultiAccess = (1+eta)*PAL;
else
  error('Unrecognized Pflag value');
end
N0effL = N0L + (2/3)*PbarMultiAccess*s.Tc*s.M;
PfullL = PIAL + N0effL*s.WFEbeta;
o.Pfull = 10*log10(PfullL);
beta2dB = s.Pbeta - o.Pfull;
beta = sqrt(10^(beta2dB/10));

%----- Simulate P
% attnFactor is the factor by which the measured power in rA and rI will be
% reduced due to s.WFE being narrower than s.WFEbeta.  It assumes a sinc^2 PSD
% profile for the desired signal.
xVec = [0:0.001:(s.WFEbeta*s.Tc/2)]';
psdVec = sinc(xVec).^2;
iidum = find(xVec <= (s.WFE*s.Tc/2));
attnFactor = sum(psdVec(iidum))/sum(psdVec);
PL = PIAL*attnFactor + N0effL*s.WFE;
o.P = 10*log10(PL) + s.sigmaP*randn(1,1);

% JJJ: I'm not sure what this is for; you didn't modify the input
% documentation to describe it.  Since my top-level script topSimPincer
% doesn't populate s.PNom, these lines cause errors when I run topSimPincer.

% PLinNoisey= 10^(o.P/10);
% PNomLin=10^(s.PNom/10);
% o.P_WRT_Nom=10*log10(PLinNoisey/PNomLin);


%----- Calculate tauHat, thetaHat
[tauHat,thetaHat] = observationError(s.etadB,s.Delta_tau,s.Delta_theta);
if(0)
figure(100);clf;
subplot(211)
plot(tauVec,real(rAIn)); grid on;
subplot(212);
plot(tauVec,imag(rAIn)); grid on;
pause
end

%----- Calculate noise in xiN
sigmaN0 = sqrt((N0L + (2/3)*PAL*s.Tc*s.M)/2/s.Ta);
% % In the case of jamming, which is indicated by s.Delta_tau >= 2, we
% % increase the effective noise floor by the amount contributed by the jamming
% % signal rI itself, treating it as another multi-access signal.  Note that for
% % large s.M (e.g., s.M >= 7), this step changes N0effL only slightly.
% if(s.Delta_tau >= 2)
%   N0effL = (2/3)*eta*PAL*s.Tc + N0effL;
% end
sigmaN = sqrt(N0effL/2/s.Ta);
if(s.taud < 0)
  error('s.taud must be non-negative');
end
if(s.taud > 1)
  Cvp = 0;
else
  Cvp = (sigmaN^2)*(1 - 2*s.taud);
end
Cvv = sigmaN^2;
C = [Cvv,Cvp;Cvp,Cvv];
Ra = chol(C);
xiNVec = Ra'*randn(2,1) + j*Ra'*randn(2,1);

%----- Calculate d
Delta_tauA = 0 - tauHat;
Delta_tauI = s.Delta_tau - tauHat;
Delta_thetaA = 0 - thetaHat;
Delta_thetaI = s.Delta_theta - thetaHat;
xiAe = sqrt(PAL)*Rcorr(-Delta_tauA - s.taud)*exp(j*Delta_thetaA);
xiAl = sqrt(PAL)*Rcorr(-Delta_tauA + s.taud)*exp(j*Delta_thetaA);
xiIe = sqrt(eta*PAL)*Rcorr(-Delta_tauI - s.taud)*exp(j*Delta_thetaI);
xiIl = sqrt(eta*PAL)*Rcorr(-Delta_tauI + s.taud)*exp(j*Delta_thetaI);
xie = beta*(xiAe + xiIe + xiNVec(1));
xil = beta*(xiAl + xiIl + xiNVec(2));
o.d = (xie - xil)/sigmaN0;

%----- Check assumption that (beta*sigmaN/sigmaN0)^2 is approximately unity
% If beta is set such that beta = 1 under the interference-free case, then
% this ratio should be near unity.
ratioCheck = (beta*sigmaN/sigmaN0)^2;



