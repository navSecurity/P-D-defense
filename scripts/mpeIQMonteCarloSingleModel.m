%% Top-level Monte-Carlo Script for simulating PD-ML Observables
% Author: Jason Gross with slight mods by Cagri Kilic

clear;clc;close all;
% set some simulation defaults
s.PA = -156;
s.etadB = 20;
s.N0 = -204;
s.Ms = 7;
s.eml = 0.3;
s.Ta = 0.1;
s.WFE = 2e6;
s.WFEbeta = 20e6;
s.Pbeta = -130.8815;
s.Tc = 1e-6;
s.sigmaP = 0.40;
s.taud = 0.15;
s.trueDelaySigma=0.03;
s.PNom=-140.0906; % parameter added to model
PNomLin=10^(s.PNom/10);
p.N = 1;
p.Tc = s.Tc;
% multipath estimator struct
%------------------------------------------------
s.nTaps=11; %select nTaps number as 5,7,15 or 41
%------------------------------------------------
mpeIn1.nTaps=s.nTaps;
mpeIn1.taps=-(s.nTaps-1)/2:(s.nTaps-1)/2;
mpeIn1.tapSize= 2/(mpeIn1.nTaps-1); % chips
mpeIn1.delayCombinations  = mpeIn1.taps' ;(-20:.01:20)';%delayGridCombinations( mpeIn.taps );
mpeIn1.invQ=inv(toeplitzQMatrix(mpeIn1.taps,mpeIn1.tapSize));



% precompute H mats
for i=1:length(mpeIn1.delayCombinations)
    delays=mpeIn1.delayCombinations(i,:);
    H=defineObservationMat(delays, mpeIn1.taps, mpeIn1.tapSize);
    mpeIn1.L(:,:,i)=pinv(H'*mpeIn1.invQ*H)*H'*mpeIn1.invQ;
    mpeIn1.Hs(:,:,i)=H;
end


Np=100000;% number of trials of theta
Nm=20;  % for each theta number of measurement vectors simulated
% store the theta for each set of measurements
thetaKeys=zeros(1,Np); %clean->0, multipath->1, spoof->2, jam->3
% assumed prior distribution
pPriorClean=.60;
pPriorMultipath=.20;
pPriorSpoof=.05;
pPriorJam=.15;

verbose=1;
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
        thetaKeys(i)=0;
        p.i = 0;
        op = simulatePincerParameters(p);
    elseif testU < pPriorClean+pPriorMultipath
        scenario = 'MULTIPATH';
        thetaKeys(i)=1;
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
        thetaKeys(i)=2;
    else
        scenario = 'JAMMING';
        p.i = 3;
        % ranging from weak jamming to severe jamming
        p.riceSdB = rand*12+3;
        p.riceSigmadB = 1;
        op = simulatePincerParameters(p);
        thetaKeys(i)=3;
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
        deltaTheta=rand*2*pi;
        s.etadB = op.etadB;
        s.Delta_tau = op.Delta_tau;
        s.Delta_theta = op.Delta_theta;
        s.Pflag = 1;
        s.colorSpec ='b.';
           
    end
    simConfigs(i)=s;
    pincerParamsOut(i)=op;
    for j=1:Nm
        
        % save 2D arrays of the sim config and PM-ML sim output
        % i is the index of the theta selection (1 to Np)
        % j is the index of the simulated measurement for each i (1 to Nm)
        [Out(i,j)] = simulateIQTaps(s);
        mpeIn1.signalTaps=Out(i,j).signalTaps;

        
        [mpeOut1(i,j)] =  multiTapMultipathEstimatorSingleModelBisect( mpeIn1);
            
        
    end
    if(verbose==1)
        waitbar(i/Np);
    end
    
end
close(h);
clear h;

if (verbose==1)
    h=waitbar(0,'Plot Progress');
end
% figure
MarkerSize=1;

cd1=[];
mp1=[];
sp1=[];
jm1=[];

cd2=[];
mp2=[];
sp2=[];
jm2=[];


cd1X=[];
mp1X=[];
sp1X=[];
jm1X=[];


cdP=[];
mpP=[];
spP=[];
jmP=[];


for i=1:Np
    
    if thetaKeys(i)==0
        
        for j=1:Nm
            PLinNoisey= 10^(Out(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            cd1=[cd1 mpeOut1(i,j).D_RSS];
            cd2=[cd2 mpeOut1(i,j).sigRSS];
            cdP=[cdP P_WRT_Nom];
            cd1X=[cd1X mpeOut1(i,j).chiSqr];

        end
        
    elseif thetaKeys(i)==1
        
        for j=1:Nm
            PLinNoisey= 10^(Out(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            mp1=[mp1 mpeOut1(i,j).D_RSS];
            mp2=[mp2 mpeOut1(i,j).sigRSS];
            mpP=[mpP P_WRT_Nom];
            mp1X=[mp1X mpeOut1(i,j).chiSqr];

        end
    elseif thetaKeys(i)==2
        
        for j=1:Nm
            PLinNoisey= 10^(Out(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            sp1=[sp1 mpeOut1(i,j).D_RSS];
            sp2=[sp2 mpeOut1(i,j).sigRSS];
            spP=[spP P_WRT_Nom];
            sp1X=[sp1X mpeOut1(i,j).chiSqr];
         
        end
    else
        
        for j=1:Nm
            PLinNoisey= 10^(Out(i,j).P/10);
            P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
            jm1=[jm1 mpeOut1(i,j).D_RSS];
            jm2=[jm2 mpeOut1(i,j).sigRSS];
            jmP=[jmP P_WRT_Nom];
            jm1X=[jm1X mpeOut1(i,j).chiSqr];
         
        end
    end
    
    waitbar(i/Np);
end


% clear h
% delete h
% save -v7.3 mpeIQMonteCarloall.mat

figure
semilogx(sp1X,spP,'r.','MarkerSize', 5);
hold on
semilogx(jm1X,jmP,'b.','MarkerSize', 5);
semilogx(mp1X,mpP,'k.','MarkerSize', 5);
semilogx(cd1X,cdP,'g.','MarkerSize', 5);
set(gca,'FontSize',16)
xlabel('RSS($$\frac{1}{a_0} exp(-j\theta_0)(\xi_k$$ - $$\hat{\xi_0}$$))','Interpreter','Latex')
ylabel('Power (dB)')
shg
axis tight
shg
title('')

 
