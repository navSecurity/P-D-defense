function [R] = Rcorr(tauVec)
% Rcorr : Return a triangular correlation function for tauVec, an N-by-1 lag
%         vector given in chips.  The standard triangular correlation function
%         is a good model for binary random sequences, such as the GPS Gold
%         codes approximate.

v = 1 - abs(tauVec(:));
M = [v';zeros(1,length(v))];
R = max(M,[],1);