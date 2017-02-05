% decisionForMinimumBayesianRiskBoundary.m create the
% P-D defense deicison regions 
% Author: Jason Gross

% note: this script assumes that the Matlab workspace
% is populated with modelPincerDistributionsMonteCarlo.m

% configuration variables
costType=0; 		% variable cost in theta 2, fixed 1, misclass 0
SAVEFIGS=1; 		% switch to save time history of detector optimzation
tau_eml=0.3; 		% spacing assumed (in chips) for the symmetric difference
power=-1.2:.4:20;	% range and resolution of power for which regions are drawn
deml=0:.4:40;		% range and resolution of symmetric differences for which regions are drawn

[demlR,powerR]=meshgrid(deml,power);

% assume a reasonable partition based on visual
% inspection of the scatter plot for the initial deicison regions

% H0: D=>[0,4) P=>[-1,1)
% H1: D=>[4,11] P=>[-1,1)
% H3:  D=>[0,4) P=>[1,20]
% H2: D=>[4,35] P=>[1,20] && D=>[11,35] P=>[-1,1]

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
        elseif (power(i)>=1 && deml(j)>=4)
            pincerDecision(i,j)=3;
        else
            pincerDecision(i,j)=4;
        end
    end
end

% plot the initial decision regoins
mymap=[0 1 0; 0 0 0; 1 0 0; 0 0 1];
figure()
colormap(mymap);
h=pcolor(deml,pow,pincerDecision);
set(h, 'EdgeColor', 'none');
xlabel('Distortion |D(eml)|')
ylabel('Power (dB)')
xlim([.05 30])
ylim([-1 19.8])
set(gca,'FontSize',14)
title('Initial Pincer Decision Partions')


%  go through all the Monte-Carlo simulated observations and determine the grid
%  location, or GridKey.  Keeping track of these will speed up the cost
%  optimzation;

%initialize variables
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
        % storing this way makes keying the data easier
        thetas(i,j)=thetaKeys(i);
        deltaThetas(i,j)=obsSimConfigs(i).Delta_theta;
        deltaTaus(i,j)=obsSimConfigs(i).Delta_tau;
        etadBs(i,j)=obsSimConfigs(i).etadB;
        
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
lastPincerDecision=pincerDecision;
% unless it is the first optimzation iteration, quit optimizing if no
% changes were made in the last iteration.
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
            end
            
	    % keep track if this location has already been evaluated
            keysAlreadyChecked=[keysAlreadyChecked key];
            
            % evaluate cost of current decision and the cost of all the other
	    % possible decisions for this  grid key location
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
                
                % 2D grid location from key value
                powLoc=mod(key,length(deml));
                demlLoc=floor((key)/length(deml));
                
                
                disp([demlLoc powLoc])
                pincerDecision(powLoc,demlLoc)=newDecision(1);
                changedDecisionIndex=changedDecisionIndex+1;

                % save the cost improvement history (for plotting)
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




