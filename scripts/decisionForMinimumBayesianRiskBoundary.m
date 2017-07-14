% create the evalution grid & evalate p(zk)incer
% Author: Jason Gross
%clc; clear; close all;
%load pincerMonteCarloData11.mat
costType=1; % variable cost in theta 2, fixed 1, misclass 0
SAVEFIGS=0;
tau_eml=0.3;
power=-1.2:.35:20;
deml=0:.35:40;
[demlR,powerR]=meshgrid(deml,power);

% start a reasonable partition based on visual
% inspection of the scatter plot for initial deicison

% clean D=>[0,4) P=>[-1,1)
% multipath D=>[4,11] P=>[-1,1)
% jam  D=>[0,4) P=>[1,20]
% spoof D=>[4,35] P=>[1,20] && D=>[11,35] P=>[-1,1]

pincerDecision=zeros(size(demlR));
for i=1:length(power)
    for j=1:length(deml)
        if(power(i)<1 && deml(j)<2)
            pincerDecision(i,j)=1;
        elseif (power(i)<1)
            if(deml(j)>=2 && deml(j)<15)
                pincerDecision(i,j)=2;
            else
                pincerDecision(i,j)=2;
            end
        elseif (power(i)>=1 && power(i)<=5 && deml(j)>=6)
            pincerDecision(i,j)=3;
        elseif (power(i)>=5 && power(i)<=10 && deml(j)>=7)
            pincerDecision(i,j)=3;
        elseif (power(i)>=10 && deml(j)>=7.5)
            pincerDecision(i,j)=3;
        else
            pincerDecision(i,j)=4;
        end
    end
end

mymap=[0 1 0; 0 0 0; 1 0 0; 0 0 1];

figure

colormap(mymap);
h=pcolor(deml,power,pincerDecision);
set(h, 'EdgeColor', 'none');
xlabel('Distortion |D(eml)|')
ylabel('Power (dB)')
xlim([.05 30])
ylim([-1 19.8])
set(gca,'FontSize',14)
title('Initial Pincer Decision Partions')


%  go through all the simulated observations and determine the grid
%  location, or GridKey.  Keeping track of these will speed up the cost
%  optimzation
dataKeys=zeros(Np,Nm);
powerSimObs=zeros(Np,Nm);
demlSimObs=zeros(Np,Nm);
thetas=zeros(Np,Nm);

deltaThetas=zeros(Np,Nm);
deltaTaus=zeros(Np,Nm);
etadBs=zeros(Np,Nm);

for i=1:Np
    for j=1:Nm
        PLinNoisey= 10^(pincerObsOut(i,j).P/10);
        P_WRT_Nom=10*log10(PLinNoisey/PNomLin);
        p=P_WRT_Nom;
        d=norm(pincerObsOut(i,j).d);
        
        [c powBinLocation] = min(abs(power-p));
        [c demlBinLocation] = min(abs(deml-d));
        
        % determine unique 1D key within the 2D grid
        dataKeys(i,j)=(demlBinLocation)*length(deml)+powBinLocation;
        powerSimObs(i,j)=p;
        demlSimObs(i,j)=d;
        % although there is only one unique theta per set of Nm obs
        % storing this way makes keying the data easier in the next step
        thetas(i,j)=thetaKeys(i);
        deltaThetas(i,j)=obsSimConfigs(i).Delta_theta;
        deltaTaus(i,j)=obsSimConfigs(i).Delta_tau;
        etadBs(i,j)=obsSimConfigs(i).etadB;
        
    end
end

% go through all sim obs & determine the cost of the sim obs and all other sim obs of the same key
% also check cost if other hypotheses were accepted
% would the cost be reduced by accepting a different hypothesis. If so, change the
% deicison.  for efficiency, keep track of the keys that have been checked

keysAlreadyChecked=[];
changedDecisionIndex=1;
if(SAVEFIGS)
    figSaveStr=sprintf('pincerDec%08d.png',changedDecisionIndex);
    saveas(gcf,figSaveStr);
end
costChange=0;
firstIteration=1;
lastPincerDecision=pincerDecision;
while(sum(sum(pincerDecision-lastPincerDecision))~=0 || firstIteration==1 )
    firstIteration=0;
    lastPincerDecision=pincerDecision;
    for i=1:Np
        for j=1:Nm
            key=dataKeys(i,j);
            % only check each grid key once
            if find(key==keysAlreadyChecked)
                continue
            end
            % determine if on a decision boundary, if not, continue
            
            onDecisionBoundary=0;
            powLoc=mod(key,length(deml));
            
            demlLoc=floor((key)/length(deml));
            if((demlLoc>1))&& (demlLoc<length(deml)) &&( (powLoc>1)&&(powLoc<length(power)))
                % disp('Here')
                if powLoc<length(power)
                    if(pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc+1,demlLoc))
                        onDecisionBoundary=1;
                    end
                end
                if powLoc<length(power) && demlLoc<length(deml)
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc+1,demlLoc+1))
                        onDecisionBoundary=1;
                    end
                end
                if(demlLoc<length(deml))
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc,demlLoc+1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc>1)
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc-1,demlLoc))
                        onDecisionBoundary=1;
                    end
                end
                if(demlLoc>1)
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc,demlLoc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(demlLoc>1 && powLoc>1)
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc-1,demlLoc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc<length(power) && demlLoc>1)
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc+1,demlLoc-1))
                        onDecisionBoundary=1;
                    end
                end
                if(powLoc>1 && demlLoc<length(deml))
                    if (pincerDecision(powLoc,demlLoc)~=pincerDecision(powLoc-1,demlLoc+1))
                        onDecisionBoundary=1;
                    end
                end
            end
            if(onDecisionBoundary==0)
                continue
            else
                disp('Evaluating')
            end
            
            keysAlreadyChecked=[keysAlreadyChecked key];
            
            % evaluate cost of current decision and all the others
            % for this grid key
            decisionOfGridKey=interp2(demlR,powerR,pincerDecision,...
                demlSimObs(i,j),powerSimObs(i,j),'nearest');
            
            thetasInGridKey=thetas(dataKeys==key);
            etadBsinGridKey=etadBs(dataKeys==key);
            deltaTauIAsInGridKey=deltaTaus(dataKeys==key);
            deltaThetaIAsInGridKey=deltaThetas(dataKeys==key);
            
            currentCost=0;
            otherDecisions=find([1,2,3,4]~=decisionOfGridKey);
            otherDecisionCosts=[0,0,0];
            
            for k=1:length(thetasInGridKey)
                %% this constant cost mat is what will need to be made variable
                currentCost=currentCost+pincerCostFunction(decisionOfGridKey,...
                    thetasInGridKey(k)+1, etadBsinGridKey(k) ,deltaTauIAsInGridKey(k), ...
                    deltaThetaIAsInGridKey(k),tau_eml,costType);
                
                for l=1:length(otherDecisions)
                    otherDecisionCosts(l)=otherDecisionCosts(l)+pincerCostFunction(otherDecisions(l),...
                        thetasInGridKey(k)+1,etadBsinGridKey(k) ,...
                        deltaTauIAsInGridKey(k), deltaThetaIAsInGridKey(k),tau_eml,costType);
                end
            end
            
            % determine if cost can be reduced, if so, accept lowest cost
            if(sum(otherDecisionCosts<currentCost))
                newDecision=otherDecisions(otherDecisionCosts==min(otherDecisionCosts));
                disp([j,i])
                % get grid location from key
                
                powLoc=mod(key,length(deml));
                %powLoc=powBinLoc+1;
                demlLoc=floor((key)/length(deml));
                
                %disp('Loc')
                %disp([demlLoc powLoc])
                pincerDecision(powLoc,demlLoc)=newDecision(1);
                changedDecisionIndex=changedDecisionIndex+1;
                %save the cost improvement history for plotting
                costChange(changedDecisionIndex)=costChange(changedDecisionIndex-1)+...
                    (currentCost-min(otherDecisionCosts));
                
                colormap(mymap);
                h=pcolor(deml,power,pincerDecision);
                set(h, 'EdgeColor', 'none');
                xlabel('Distortion |D(eml)|')
                ylabel('Power (dB)')
                xlim([.05 30])
                set(gca,'FontSize',14)
                if(SAVEFIGS)
                    figSaveStr=sprintf('pincerDec%08d.png',changedDecisionIndex);
                    saveas(gcf,figSaveStr);
                end
            end
        end
    end
end




