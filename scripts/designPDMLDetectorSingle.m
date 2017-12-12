% decisionForMinimumBayesianRiskBoundary.m create the
% PDML defense decision regions
% Author: Jason Gross with slight mods by Cagri Kilic
% -------------------------------------------------------------------------
% NOTE: this script assumes that the Matlab workspace
% is populated with mpeIQMonteCarlo.m if not please uncomment below!
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% clear all;
% load mpeIQMonteCarloall.mat
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% configuration variables
% -------------------------------------------------------------------------
nTaps=11;            % select the taps number as 5,7,15 or 41 if you have
% -------------------------------------------------------------------------
% variable cost in theta 2, fixed 1, misclass 0
SAVEFIGS=0;             % switch to save time history of detector optimzation
tau_eml=2/(nTaps-1);    % spacing assumed (in chips) for the symmetric difference
pow=-1.5:.3:20;        % range and resolution of power for which regions are drawn
d1=[01:5:100 101:50:500 550:500:8500];		% range and resolution of symmetric differences for which regions are drawn
[da1R,powerR]=meshgrid(d1,pow);

% assume a reasonable partition based on visual
% inspection of the scatter plot for the initial decision regions

pdmlDecision=zeros(size(da1R));
benignSpoofMat=zeros(size(da1R));

%% 41taps
for i=1:length(pow)
    for j=1:length(d1)
        if (pow(i)<1.1&& d1(j)<30  )
            pdmlDecision(i,j)=1;
        elseif (pow(i)<1.1 && d1(j)>=30 )
            pdmlDecision(i,j)=2;
        elseif (pow(i)>=1.1&& d1(j)<50)
            pdmlDecision(i,j)=4;
        else
            pdmlDecision(i,j)=3;
        end
    end
end
% plot the initial decision regoins
mymap=[0 1 0; 0 0 0; 1 0 0; 0 0 1];

figure()
colormap(mymap);
h=pcolor(log10(d1),pow,pdmlDecision);
set(h, 'EdgeColor', 'none');
%xlabel('Distortion |a1|')
%ylabel('Power (dB)')
set(gca,'FontSize',14)
title('Initial Pincer Decision Partions')
drawnow
%saveas(gcf,'pincerDecinit.png');
%saveas(gcf,'pincerDecinit.eps');
%for a=0:2
a=1;
costType=a;

%initialize variables
dataKeys=zeros(Np,Nm);
powerSimObs=zeros(Np,Nm);
da0a1SimObs=zeros(Np,Nm);
thetas=zeros(Np,Nm);

deltaThetas=zeros(Np,Nm);
deltaTaus=zeros(Np,Nm);
etadBs=zeros(Np,Nm);
for i=1:Np
    for j=1:Nm
        PLinNoisey= 10^(Out(i,j).P/10);
        P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
        p=P_WRT_Nom;
        d=(mpeOut1(i,j).chiSqr);
        
        [c powBinLocation] = min(abs(pow-p));
        [c da0a1BinLocation] = min(abs(d1-d));
        
        % determine unique 1D key within the 2D grid
        dataKeys(i,j)=(da0a1BinLocation)*length(d1)+powBinLocation;
        powerSimObs(i,j)=p;
        da0a1SimObs(i,j)=d;
        % although there is only one unique theta per set of Nm obs
        % storing this way makes keying the data easier
        thetas(i,j)=thetaKeys(i);
        deltaThetas(i,j)=simConfigs(i).Delta_theta;
        deltaTaus(i,j)=simConfigs(i).Delta_tau;
        etadBs(i,j)=simConfigs(i).etadB;
        
    end
end

% go through all the simulated obs & determine if it lies on a decision
% region boundary. if so, evaluate the current cost
% and find all other simulated obs of the same key location.
% determine if the cost of the location would be lower
% if another hypotheses were accepted.
% if so, change the for the location decision accordingly.
% for efficiency, keep track of the keys that have been checked and
% only check each key them once.
% iterate until now more boundaries are changed.

keysAlreadyChecked=[];
changedDecisionIndex=1;
if(SAVEFIGS)
    figSaveStr=sprintf('pincerDec%08d.png',changedDecisionIndex);
    saveas(gcf,figSaveStr);
end

costChange=0;
firstIteration=1;
lastPincerDecision=pdmlDecision;
% unless it is the first optimization iteration, quit optimizing if no
% changes were made in the last iteration.
while(sum(sum(pdmlDecision-lastPincerDecision))~=0 || firstIteration==1 )
    firstIteration=0;
    lastPincerDecision=pdmlDecision;
    for i=1:Np
        for j=1:Nm
            key=dataKeys(i,j);
            % only check each grid key once
            if find(key==keysAlreadyChecked)
                continue
            end
            % determine if on a decision boundary, if not, continue
            onDecisionBoundary=0;
            powLoc=mod(key,length(d1));
            
            da1Loc=floor((key)/length(d1));
            if((da1Loc>1))&& (da1Loc<length(d1)) &&( (powLoc>1)&&(powLoc<length(pow)))
                
                if powLoc<length(pow)
                    if(pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc+1,da1Loc))
                        onDecisionBoundary=1;
                    end
                end
                if powLoc<length(pow) && da1Loc<length(d1)
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc+1,da1Loc+1))
                        onDecisionBoundary=1;
                    end
                end
                if(da1Loc<length(d1))
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc,da1Loc+1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc>1)
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc-1,da1Loc))
                        onDecisionBoundary=1;
                    end
                end
                if(da1Loc>1)
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc,da1Loc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(da1Loc>1 && powLoc>1)
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc-1,da1Loc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc<length(pow) && da1Loc>1)
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc+1,da1Loc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc>1 && da1Loc<length(d1))
                    if (pdmlDecision(powLoc,da1Loc)~=pdmlDecision(powLoc-1,da1Loc+1))
                        onDecisionBoundary=1;
                    end
                end
            end
            
            if(onDecisionBoundary==0)
                continue
            end
            
            % keep track if this location has already been evaluated
            keysAlreadyChecked=[keysAlreadyChecked key];
            
            % evaluate cost of current decision and the cost of all the other
            % possible decisions for this  grid key location
            decisionOfGridKey=interp2(da1R,powerR,pdmlDecision,...
                da0a1SimObs(i,j),powerSimObs(i,j),'nearest');
            
            thetasInGridKey=thetas(dataKeys==key);
            etadBsinGridKey=etadBs(dataKeys==key);
            deltaTauIAsInGridKey=deltaTaus(dataKeys==key);
            deltaThetaIAsInGridKey=deltaThetas(dataKeys==key);
            
            currentCost=0;
            otherDecisions=find([1,2,3,4]~=decisionOfGridKey);
            otherDecisionCosts=[0,0,0];
            powLoc=mod(key,length(d1));
            da1Loc=floor((key)/length(d1));
            
            
            for k=1:length(thetasInGridKey)
                
                currentCost=currentCost+ pdmlCostFunctionSingle(decisionOfGridKey,...
                    thetasInGridKey(k)+1, etadBsinGridKey(k) ,deltaTauIAsInGridKey(k), ...
                    deltaThetaIAsInGridKey(k),tau_eml,costType);
                %                  benignSpoofMat(powLoc,da1Loc)=benignSpoof(decisionOfGridKey,...
                %                     thetasInGridKey(k)+1, etadBsinGridKey(k) ,deltaTauIAsInGridKey(k), ...
                %                     deltaThetaIAsInGridKey(k),tau_eml,costType);
                
                for l=1:length(otherDecisions)
                    otherDecisionCosts(l)=otherDecisionCosts(l)+ pdmlCostFunctionSingle(otherDecisions(l),...
                        thetasInGridKey(k)+1,etadBsinGridKey(k) ,...
                        deltaTauIAsInGridKey(k), deltaThetaIAsInGridKey(k),tau_eml,costType);
                end
            end
            
            % determine if cost can be reduced, if so, accept lowest cost
            if(sum(otherDecisionCosts<currentCost))
                newDecision=otherDecisions(otherDecisionCosts==min(otherDecisionCosts));
                
                % 2D grid location from key value
                powLoc=mod(key,length(d1));
                da1Loc=floor((key)/length(d1));
                
                
                disp([da1Loc powLoc])
                pdmlDecision(powLoc,da1Loc)=newDecision(1);
                changedDecisionIndex=changedDecisionIndex+1;
                
                % save the cost improvement history (for plotting)
                costChange(changedDecisionIndex)=costChange(changedDecisionIndex-1)+...
                    (currentCost-min(otherDecisionCosts));
            end
        end
    end
end
[M,N]=size(pdmlDecision);
for i=2:M-1
    for j=2:N-1
        
        decision=pdmlDecision(i,j);
        boundaryDecision=pdmlDecision(i,j+1);
        boundaryCount=0;
        if ((pdmlDecision(i,j-1)==boundaryDecision))
            boundaryCount= boundaryCount+1;
        end
        if ((pdmlDecision(i-1,j)==boundaryDecision))
            boundaryCount= boundaryCount+1;
        end
        if ((pdmlDecision(i+1,j)==boundaryDecision))
            boundaryCount= boundaryCount+1;
        end
        if (decision~=boundaryDecision)&&(boundaryCount>0)
            disp('Found Island')
            pdmlDecision(i,j)=boundaryDecision;
        end
    end
end

figure()
colormap(mymap);
h=pcolor(log10(d1),pow,pdmlDecision+benignSpoofMat);
set(h, 'EdgeColor', 'none');
%xlabel('Distortion |a1|')
%ylabel('Power (dB)')
set(gca,'FontSize',14)

figSaveStrpng=sprintf('pincerDec%d.png',a);
saveas(gcf,figSaveStrpng);
figSaveStrfig=sprintf('pincerDec%d.fig',a);
saveas(gcf,figSaveStrfig);
figSaveStrfig=sprintf('pincerDec%d.eps',a);
saveas(gcf,figSaveStrfig);

if costType==0
    save detectionCF0a1.mat d1 pow pdmlDecision
elseif costType==1
    save detectionCF1a1.mat d1 pow pdmlDecision
elseif costType==2
    save detectionCF2a1.mat d1 pow pdmlDecision
end
%end
