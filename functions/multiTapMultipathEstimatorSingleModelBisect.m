function [ output ] = multiTapMultipathEstimatorSingleModelBisect( input )
% Author: Jason N. Gross
for i=1:length(input.delayCombinations)
    X(:,i)=input.L(:,:,i)*input.signalTaps;
    dz=(input.signalTaps-input.Hs(:,:,i)*X(:,i));
    err(i)=norm(dz'*input.invQ*dz);
end

[minErr,minInd1]=min(err);
secondSmallestError=min(err(err~=minErr));
minInd2=find(err==secondSmallestError);
firstIter=1;
bisectErr=minErr;
bestDelay=input.delayCombinations(minInd1);
secondBestDelay=input.delayCombinations(minInd2);
iterCount=0;
lastXX=X(:,minInd1);
while(((((minErr-bisectErr)>1e-30) )|| firstIter )&&(iterCount<100))
    
    if(~firstIter)
        if(abs(bisectDelay-bestDelay)<=abs(bisectDelay-secondBestDelay))
            secondBestDelay=bestDelay;
        end
        bestDelay=bisectDelay;
        minErr=bisectErr;
        lastXX=XX;
    end
    firstIter=0;
    bisectDelay=(bestDelay+secondBestDelay)/2;
    
    H=defineObservationMat(bisectDelay, input.taps, input.tapSize);
    
    
    L=pinv(H'*input.invQ*H)*H'*input.invQ;
    
    XX=L*input.signalTaps;
    dz=(input.signalTaps-H*XX);
    bisectErr=norm(dz'*input.invQ*dz);
    if(bisectErr>minErr && bisectErr < secondSmallestError)
     secondSmallestError=bisectErr;
     secondBestDelay=bisectDelay;
     firstIter=1;
 
   
    end

    iterCount=iterCount+1;
end

output.chiSqr=bisectErr;
output.estDelay=bisectDelay;
output.a0=norm(XX);
output.theta0=atan2(imag(XX),real(XX));

for i=1:length(input.taps)
    output.sigEst(i,1)=output.a0*exp(1i*output.theta0)*GPSCAAutoCorr((input.taps(i)-output.estDelay(1))*input.tapSize);
   
end
%
output.Rhat=(1/output.a0)*exp(-1i*output.theta0)*(input.signalTaps-output.sigEst(:,1));
output.sigHat=exp(-1i*output.theta0)*(input.signalTaps-output.sigEst(:,1));

output.D_RSS=output.Rhat'*output.Rhat;
output.sigRSS=output.sigHat'*output.sigHat;



