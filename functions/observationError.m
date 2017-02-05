function [tauError,thetaError] = observationError(etadB,Delta_tau,Delta_theta)
% observationError : Return code- and carrier-phase measurement error (in
%                    chips and radians, respectively), as caused by the
%                    interference with power advantage etadB (in dB), code
%                    offset Delta_tau (in chips), and phase offset Delta_theta
%                    (in radians).
%

eta = 10^(etadB/10);
tau_stride = 0.001;  % in chips
tauVec = [(-1 - abs(Delta_tau)):tau_stride:(1 + abs(Delta_tau))]';
rAn = Rcorr(tauVec);
rIn = sqrt(eta)*Rcorr(tauVec - Delta_tau)*exp(j*Delta_theta);
rAIn = rAn + rIn;
[mx,iidum] = max(abs(rAIn));
tauError = abs(tauVec(iidum));
thetaError = angle(rAIn(iidum));
