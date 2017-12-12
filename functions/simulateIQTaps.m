function [o] = simulateIQTaps(s)
% 			   : Simulate the multiple taps observables used for
%                             hypothesis testing in the PD-ML defense.
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
%      Delta_tau = tau_I - tau_A, the difference between the interference
%                    and the authentic code phases (chips)
%
%    Delta_theta = theta_I - theta_A, the difference between the
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
%               TC = Chip interval (seconds)
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
%   trueDelaySigma standard deviation of autehtical delay (chips)
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
%     beta = AGC scaling factor. The incoming signal r is scaled
%            by beta such that the power in r*beta is Pbeta (dBW). In the
%            interference-free case, beta = 1.
%
%   sigmaN = Standard deviation of correlation function due to thermal and
%            multi-access noise plus interference
%
%  sigmaN0 = The value of sigmaN for the interference-free case
%

%+------------------------------------------------------------------------------+
% References: 
%
% Author:  This is based on Todd Humphreys simulatedPincerObservables and
%          modified by Jason Gross to simulated n-Taps
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
N0effL = N0L + (2/3)*PbarMultiAccess*s.Tc*s.Ms;
PfullL = PIAL + N0effL*s.WFEbeta;
o.Pfull = 10*log10(PfullL);

% s.Pbeta should be set such that the interference-free case yields beta = 1
beta2dB = s.Pbeta - o.Pfull;
o.beta = sqrt(10^(beta2dB/10));
ntaps=s.nTaps;
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
o.CN0=(10*log10(PL))-s.N0;
PLinNoisey= 10^(o.P/10);
PNomLin=10^(s.PNom/10);
o.P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
%----- Calculate tauHat, thetaHat
delta_tau = 2/(ntaps-1);
tauVec = [-1:delta_tau:1]';
%% authentic phase and delay
theta=2*pi*rand;
delay=s.trueDelaySigma*randn; % assume some authentic tracking delay
rAn = sqrt(PAL)*Rcorr(tauVec+delay)*exp(1i*theta);
rIn = sqrt(eta*PAL)*Rcorr(tauVec +delay- s.Delta_tau)*exp(1i*(theta-s.Delta_theta));
rAIn = rAn + rIn;
[mx,iidum] = max(abs(rAIn));
tauHat = tauVec(iidum);
thetaHat = angle(rAIn(iidum));
if(0)
figure(100);clf;
subplot(211)
plot(tauVec,real(rAIn),'-o'); grid on;
subplot(212);
plot(tauVec,imag(rAIn),'-o'); grid on;
disp(length(tauVec))
%pause
end

%----- Calculate noise in xiN
o.sigmaN0 = sqrt((N0L + (2/3)*PAL*s.Tc*s.Ms)/2/s.Ta);
% % In the case of jamming, which is indicated by s.Delta_tauIA >= 2, we
% % increase the effective noise floor by the amount contributed by the jamming
% % signal rI itself, treating it as another multi-access signal.  Note that for
% % large s.M (e.g., s.M >= 7), this step changes N0effL only slightly.
% if(s.Delta_tauIA >= 2)
%   N0effL = (2/3)*eta*PAL*s.Tc + N0effL;
% end
o.sigmaN = sqrt(N0effL/2/s.Ta);
if(s.taud < 0)
  error('s.taud must be non-negative');
end
if(s.taud > 1)
  Cvp = 0;
else
  Cvp = (o.sigmaN^2)*(1 - 2*s.taud);
end
Cvv = o.sigmaN^2;
C = [Cvv,Cvp;Cvp,Cvv];
Ra = chol(C);
xiNVec = Ra'*randn(2,1) + 1i*Ra'*randn(2,1);
CC=zeros(ntaps,ntaps);
for jj=1:ntaps
    for kk=1:ntaps
        if (jj==(ntaps+1)/2 && kk==(ntaps+1)/2)
            CC(jj,kk)=o.sigmaN^2;
        else
            jChips=((ntaps+1)/2-jj)/(ntaps-1);
            kChips=((ntaps+1)/2-kk)/(ntaps-1);
            CC(jj,kk)=(o.sigmaN^2)*(1 - abs(jChips-kChips));
        end
    end
end

Rtaps = chol(CC);%/sqrt(ntaps);


N=Rtaps'*randn(ntaps,1);
xiNVecTaps = Rtaps'*randn(ntaps,1) + 1i*Rtaps'*randn(ntaps,1);

o.signalTaps=o.beta*(rAIn'+xiNVecTaps)/o.sigmaN0;
if(0)
figure(200);clf;
subplot(211)
plot(tauVec,real(o.signalTaps),'-o'); grid on;
subplot(212);
plot(tauVec,imag(o.signalTaps),'-o'); grid on;

%pause
end
%----- Calculate d
Delta_tauA = 0 - tauHat;
Delta_tauI = s.Delta_tau - tauHat;
Delta_thetaA = 0 - thetaHat;
Delta_thetaI = s.Delta_theta - thetaHat;
xiAe = sqrt(PAL)*Rcorr(-Delta_tauA - s.taud)*exp(1i*Delta_thetaA);
xiAl = sqrt(PAL)*Rcorr(-Delta_tauA + s.taud)*exp(1i*Delta_thetaA);
xiIe = sqrt(eta*PAL)*Rcorr(-Delta_tauI - s.taud)*exp(1i*Delta_thetaI);
xiIl = sqrt(eta*PAL)*Rcorr(-Delta_tauI + s.taud)*exp(1i*Delta_thetaI);
xie = o.beta*(xiAe + xiIe + xiNVec(1));
xil = o.beta*(xiAl + xiIl + xiNVec(2));
o.d = (xie - xil)/o.sigmaN0;

%----- Check assumption that (beta*sigmaN/sigmaN0)^2 is approximately unity
% If beta is set such that beta = 1 under the interference-free case, then
% this ratio should be near unity.
ratioCheck = (o.beta*o.sigmaN/o.sigmaN0)^2;

% Return correlation function for tau, an N-by-1 lag vector given in chips
function [R] = Rcorr(tau)
v = 1 - abs(tau(:));
M = [v';zeros(1,length(v))];
R = max(M,[],1);



