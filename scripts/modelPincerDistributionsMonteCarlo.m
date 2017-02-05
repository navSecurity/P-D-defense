% modelPincerDistributionsMonteCarlo.m: Monte-Carlo Script for simulating pincer observables
% the output of this script is required for determining decision regions with 
% decisionForMinimumBayesianRiskBoundary.m
% Author: Jason Gross

clear;clc;close all;

% output file name
monteCarloFileName='pdMonterCarlo.mat'

% set some simulation defaults for the pincer sim
s.PA = -156;
s.N0 = -204;
s.M = 7;
s.taud = 0.15;
s.Ta = 0.1;
s.WFE = 2e6;
s.WFEbeta = 20e6;
s.Pbeta = -129.8583;
s.Tc = 1e-6;
s.sigmaP = 0.45;
PNom=-139.76; % parameter added to model
PNomLin=10^(PNom/10);
p.N = 1;
p.Tc = s.Tc;


% monte-carlo configuration variables
Np=10000;% number of trials of theta
Nm=20; % for each theta number of measurement vectors simulated
verbose=1;  %  used for waitbars

% store the theta for each set of measurements
thetaKeys=zeros(1,Np); %H0->0, H1->1, H2->2, H3->3

% prior for each hypothesis
pPriorClean=.60;
pPriorMultipath=.20;
pPriorSpoof=.05;
pPriorJam=.15;



if (verbose==1)
    h=waitbar(0,'Simulation Progress');
end

for i=1:Np
    clearvars p
    p.N = 1;
    p.Tc = s.Tc;
    testU = rand;
    
    if testU < pPriorClean
        scenario = 'CLEAN';
        p.i = 0;
        op = simulatePincerParameters(p);
    elseif testU < pPriorClean+pPriorMultipath
        scenario = 'MULTIPATH';
        p.i = 1;
        p.elDeg = raylrnd(3)+5;
        s.Pflag = 0;
        op = simulatePincerParameters(p);
    elseif testU < pPriorClean+pPriorMultipath+pPriorSpoof
        scenario = 'SPOOFING';
        p.i = 2;
        p.elDeg = 3;
        p.muEtadB =1; 
        p.sigmaEtadB =1;
        op = simulatePincerParameters(p);
    else
        scenario = 'JAMMING';
        p.i = 3;
        % ranging from weak jamming to severe jamming
        p.riceSdB = rand*12+3;
        p.riceSigmadB = 1;
        op = simulatePincerParameters(p);
        
    end
    
    
    if strcmp(scenario,'CLEAN')
        s.WFE = 2e6;
        s.WFEbeta= 10e6;
        s.Pbeta=  -131.9641;
        s.etadB = op.etadB;
        s.Delta_tau = op.Delta_tau;
        s.Delta_theta = op.Delta_theta;
        s.Pflag = 0;
        s.colorSpec ='g.';
        
        
    elseif strcmp(scenario,'MULTIPATH')
        s.WFE = 2e6;
        s.WFEbeta= 10e6;
        s.Pbeta=  -131.9641;
        s.etadB = op.etadB;
        s.Delta_tau = op.Delta_tau;
        s.Delta_theta = op.Delta_theta;
        s.Pflag = 0;
        s.colorSpec ='k.';
        
    elseif strcmp(scenario,'SPOOFING')
        s.WFE = 2e6;
        s.WFEbeta= 20e6;
        s.Pbeta=  -128.9726;
        s.etadB = op.etadB;
        s.Delta_tau = op.Delta_tau;
        s.Delta_theta = op.Delta_theta;
        s.Pflag = 2;
        s.colorSpec ='r.';
    elseif strcmp(scenario,'JAMMING')
        
        s.WFE = 2e6;
        s.WFEbeta= 10e6;
        s.Pbeta=  -131.9641;
        s.etadB = op.etadB;
        s.Delta_tau = op.Delta_tau;
        s.Delta_theta = op.Delta_theta;
        s.Pflag = 1;
        s.colorSpec ='b.';
        
        
    end
    obsSimConfigs(i)=s;
    thetaKeys(i)=p.i;
    pincerParamsOut(i)=op;
    for j=1:Nm
        
        % save 2D arrays of the sim config and pincer sim output
        % i is the index of the theta selection (1 to Np)
        % j is the index of the simulated measurement for each i (1 to Nm)
        [pincerObsOut(i,j)] = simulatePincerObservables(s);
        
    end
    if(verbose==1)
        waitbar(i/Np);
    end
    
end
close(h);
clear h;
save -v7.3 monteCarloFileName

% scatter plot
if (verbose==1)
    h=waitbar(0,'Plot Progress');
end
figure
MarkerSize=1;
cdD=[];
cdP=[];
mpD=[];
mpP=[];
spD=[];
spP=[];
jmD=[];
jmP=[];
for i=1:Np
    
    if thetaKeys(i)==0
        
        for j=1:Nm
            PLinNoisey= 10^(pincerObsOut(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            cdD=[cdD (norm(pincerObsOut(i,j).d))];
            cdP=[cdP P_WRT_Nom];
        end
        
    elseif thetaKeys(i)==1
        
        for j=1:Nm
            PLinNoisey= 10^(pincerObsOut(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            mpD=[mpD (norm(pincerObsOut(i,j).d))];
            mpP=[mpP P_WRT_Nom];
        end
    elseif thetaKeys(i)==2
        
        for j=1:Nm
            PLinNoisey= 10^(pincerObsOut(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            spD=[spD (norm(pincerObsOut(i,j).d))];
            spP=[spP P_WRT_Nom];
        end
    else

        for j=1:Nm
            PLinNoisey= 10^(pincerObsOut(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            jmD=[jmD (norm(pincerObsOut(i,j).d))];
            jmP=[jmP P_WRT_Nom];
        end
        
        
        
    end
    
    if(verbose==1)
    	waitbar(i/Np);
    end
end
figure
plot(spD,spP,'r.','MarkerSize', 1);
hold on
plot(jmD,jmP,'b.','MarkerSize', 1);
plot(mpD,mpP,'k.','MarkerSize', 1);
plot(cdD,cdP,'g.','MarkerSize', 1);
hold on

xlim([0 30])
ylim([-1 20])
xlabel('Distortion D(eml)')
set(gca,'FontSize',14)
 ylabel('Power (dB)')
 if(verbose==1)
    close(h);
end




