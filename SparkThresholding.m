function [bw,Threshold]=SparkThresholding(I,varargin)
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
    p.addParameter('ScaleN',[2,1,6],@(x)(numel(x)==3 && min(x)>0 && x(1)<=x(3)));
    parse(p, varargin{:});
    p=p.Results;

    ScaleN=p.ScaleN;
    Threshold=p.Threshold;

    %% Thresholding here.
    fprintf('%3.0f%%',0);
    StepT=ScaleN(1):ScaleN(2):ScaleN(3);
    StepTNum=numel(StepT);
    siz=size(I);
    
    bw=false(siz);
    for k=1:StepTNum
        h=[ones(StepT(k),1);-ones(StepT(k),1)];
        cA=imfilter(I,reshape(h,[1,1,numel(h)]),'symmetric','same','conv');
        variation=mad(cA,1,3)*1.4826;
        currbw=bsxfun(@ge,cA,variation.*Threshold);
        shiftLen=StepT(k)-1;      
        if shiftLen>=1
            currbw=cat(3,false(siz(1),siz(2),shiftLen),currbw(:,:,1:(siz(3)-shiftLen)));
        end
        bw=bw+currbw;
        clear('cA','currbw','shiftLen')
        fprintf('\b\b\b\b%3.0f%%',floor(k/StepTNum*100));
    end
    % bw=bw>=1;
end



% function [bw,Threshold]=SparkThresholding(I,varargin)
%     % Thresholding sparks.
%     %   Input:
%     %       ScaleN:     2-6 points with 1 as step in the third dimension.
%     %       Threshold:  default=3.
%     %   Output:
%     %       bw: the binayr image stack marking positive pixels.
%     
%     %% Inputs.
%     p=inputParser;
%     p.addParameter('Threshold',[],  @(x)(isscalar(x) && x>0) || isempty(x) ||...
%         (size(x,1)==size(I,1)&&size(x,2)==size(I,2)));
%     p.addParameter('ScaleN',[2,1,6],@(x)(numel(x)==3 && min(x)>0 && x(1)<=x(3)));
%     parse(p, varargin{:});
%     p=p.Results;
% 
%     ScaleN=p.ScaleN;
%     Threshold=p.Threshold;
%     
% 
%     %% Thresholding here.
%     fprintf('%3.0f%%',0);
%     StepT=ScaleN(1):ScaleN(2):ScaleN(3);
%     StepTNum=numel(StepT);
%     siz=size(I);
%     
%     bw=false(siz);
%     for k=1:StepTNum
%         h=[ones(StepT(k),1);-ones(StepT(k),1)];
%         cA=imfilter(I,reshape(h,[1,1,numel(h)]),'symmetric','same','conv');
%         variation=mad(cA,1,3)*1.4826;
%         % variation=QuadR(cA)*1.5;
%         % if k==1; DataCV=variation; end
%         currbw=bsxfun(@ge,cA,variation.*Threshold);
%         %Optional: currbw=cat(3,false(siz(1),siz(2),1),currbw(:,:,1:(siz(3)-1)));
%         shiftLen=StepT(k)-1; % shiftLen=floor(StepT(k)/2)-1;        
%         if shiftLen>=1
%             currbw=cat(3,false(siz(1),siz(2),shiftLen),currbw(:,:,1:(siz(3)-shiftLen)));
%         end
%         bw=bw+currbw;
%         clear('cA','currbw','shiftLen')
%         fprintf('\b\b\b\b%3.0f%%',floor(k/StepTNum*100));
%     end
%     bw=bw>=1;
% end
% function variation=QuadR(cA)
%     cA=sort(cA,3,'ascend');
%     Isiz=size(cA,3);
%     variation=cA(:,:,round(Isiz/4*3))-cA(:,:,round(Isiz/2))*2;
% end
%%
% function Threshold=AdaptiveThreshold(DataCV)
%     % TH=@(offset,Plateau,tau,x_shift,x) offset+(Plateau-offset)*(1-exp(-1/tau*(x-x_shift)));
%     p=[5.85753,0.31154,17.04231,5.86437];
%     TH=@(x,p) p(4)+(p(1)-p(4))*(1-exp(-p(2)*(x-p(3)))) + 0.3650;
%     Threshold=TH(DataCV,p);
%     
%     Threshold(isnan(Threshold))=TH(inf,p);
%     Threshold(DataCV<0)=TH(0,p);
% end
