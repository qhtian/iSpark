function [peaks]=SparkSpatialThresholding(I,varargin)
    % Thresholding sparks.
    %   Input:
    %       ScaleN:     2-6 points with 1 as step in the third dimension.
    %       Threshold:  default=3.
    %   Output:
    %       bw: the binayr image stack marking positive pixels.
    
    %% Inputs.
    p=inputParser;
    p.addParameter('Threshold',[],  @(x)(isscalar(x) && x>0) || isempty(x) ||...
        (size(x,1)==size(I,1)&&size(x,2)==size(I,2)));
    p.addParameter('xyt_dim',[0.28,0.28,10],@(x)(numel(x)==3 && min(x)>0 && x(1)<=x(3)));
    parse(p, varargin{:});
    p=p.Results;
    
    ScaleN=p.xyt_dim;
    ScaleN=round(1.12/ScaleN(1));

    
    %% Thresholding here.
    fprintf('%3.0f%%',0);
    StepT=ScaleN; %(1):ScaleN(2):ScaleN(3);
    StepTNum=numel(StepT);
    
    peaks=[];
    I=bsxfun(@minus,I,median(I,3));
    for k=1:StepTNum
        if isempty(peaks)
            peaks=calCORRX(I,StepT(k));
            peaks=peaks+calCORRY(I,StepT(k));
        else
            peaks=peaks+calCORRX(I,StepT(k));
            peaks=peaks+calCORRY(I,StepT(k));
        end
      
        fprintf('\b\b\b\b%3.0f%%',floor(k/StepTNum*100));
    end
end


function peaks=calCORRX(I,StepT)
    h=[ones(StepT,1);-ones(StepT,1)];
    cA=imfilter(I,h','symmetric','same','conv');
    peaks=zeros(size(cA));
    
    bw=cA>0;
    cApos=cA.*bw;
    cAneg=-cA.*(~bw);
    for k=1:StepT*2
        a=circshift(cApos,[0,k]);
        b=circshift(cAneg,[0,-k]);
        peaks=peaks+(a.*b);
        disp(corr(a(:),b(:)))
    end    
end

function peaks=calCORRY(I,StepT)
    h=[ones(StepT,1);-ones(StepT,1)];
    cA=imfilter(I,h,'symmetric','same','conv');
    peaks=zeros(size(cA));
    
    bw=cA>0;
    cApos=cA.*bw;
    cAneg=-cA.*(~bw);
    for k=1:StepT*2
        a=circshift(cApos,[k,0]);
        b=circshift(cAneg,[-k,0]);
        peaks=peaks+(a.*b);
    end    
end