function [CellMask,algo]=CellMasking(Data,MaskAlgorithm,xy_dim)
    % Detect cell mask, and select the largest object to analyze.
    % For this, the syntax is: [CellMask,~]=CellMasking(Data,,FineMask,LiveView);
    %
    % If the "OutOfNoise" is used to calculate the mask, give it a theshold
    % value, e.g. 2 times of noise sigma. The syntax is:
    %                          [CellMask,~]=CellMasking(Data,2,LiveView,);
    %
    % If you want to see the detected mask immediately, use a number
    % in seconds that makes LiveView>0.
    
    if nargin==1
        MaskAlgorithm=1;
    end
    
    if size(Data,3)>1
        % if MaskAlgorithm==4
        %     Data=std(single(Data),0,3);
        % else
        Data=mean(Data,3);
        % end
    end
    % dData=smoothn(Data);
    
    switch MaskAlgorithm
        case 1
            CellMask = Masking_Isodata(Data);
            algo='Isodata algorithm';
        case 2
            if nargin<=2
                disp('  --- No enough inputs for local adaptive thresholding, using Msize=9 instead. ---');
                Msize=9;
            else
                Msize=ceil(2/mean(xy_dim));
            end
            
            CellMask = neuronAxonIsolation(Data,Msize);
            algo='Local adaptive TH';
        case 3
            CellMask = Masking_Isodata(Data);
            Ibw=sort(Data(~CellMask),'ascend');
            Ibw_center=Ibw(round(numel(Ibw)*0.5));
            Ibw_sigma=Ibw_center-Ibw(round(numel(Ibw)*0.159));
            CellMask=Data>=Ibw_center+Ibw_sigma*3;
            algo='Isolated sampleBgr';
        case 4
            CellMask = true(size(Data,1),size(Data,2));
            algo='Cell Masking skipped';
    end
end


%% Isodata segmentation.
function CellMask=Masking_Isodata(I)
    if size(I,3)>1; I=mean(I,3); end
    if ~isfloat(I); I=double(I); end

    try
        Isort=sort(I(:));    Imin=Isort(1);    Imax=Isort(end);
        step=(Isort(round(numel(Isort)*0.75))-Isort(round(numel(Isort)*0.25)))/1000;
        clear('Isort');
        
        x=Imin:step:Imax;    [y,x]=hist(I(:),x);
        x=reshape(x,[numel(x),1]);    y=reshape(y,[numel(y),1]);
        
        % Apply 0~99.5% range.
        y_sum=cumsum(y);bw=y_sum<y_sum(end)*0.995;
        x=x(bw);    y=y(bw);    clear('y_sum','bw');
        
        Imax=max(x);  Mu=(Imax+Imin)/2;  Mu=double(Mu);
        options=optimset('MaxIter',100000000,'Display','off','TolX',1e-4);
        p=fminsearch(@(mu) yCenter(y,x,mu),Mu,options);
        
        % Mask here.
        CellMask=I>p;
        % CellMask=imdilate(CellMask,[0,1,0;1,1,1;0,1,0]);
        CellMask=imfill(CellMask,'holes');
        % CellMask=imerode(CellMask,[0,1,0;1,1,1;0,1,0]);
        
        % Find out the largest area to serve as a cell.
        CellMaskLabel = bwconncomp(CellMask,ones(3,3));
        area = cellfun(@numel, CellMaskLabel.PixelIdxList);
        idxToKeep=area>=max(area*0.5);
        CellMaskLabel.PixelIdxList=CellMaskLabel.PixelIdxList(idxToKeep);
        CellMaskLabel.NumObjects=sum(idxToKeep);
        
        CellMask=false(size(I));
        for k = 1 : CellMaskLabel.NumObjects
            CellMask(CellMaskLabel.PixelIdxList{k}) = true;
        end
    catch
        CellMask=false(size(I));
    end
end

function funevl=yCenter(y,x,mu)
    bw1=x<mu;
    C1=sum(y.*(x.*bw1))/sum(y.*bw1);
    bw2=~bw1;
    C2=sum(y.*(x.*bw2))/sum(y.*bw2);
    funevl=sum(((x(bw1)-C1).*y(bw1)).^2)+sum(((x(bw2)-C2).*y(bw2)).^2);
end


%% Local adaptive segmentation for neuron cells.
function bw=neuronAxonIsolation(I,Msize)
    I=single(I);
    
    % Set threshold.
    threshold=2;
    
    % Cycle to calculate adaptive local average.
    h=fspecial('disk',Msize);
    I_last=I;
    bw=false(size(I));
    bw_sum_diff=inf;
    k=0;
    while (bw_sum_diff>numel(I)*3.6808e-05) && (k<=500)
        k=k+1;
        I_mean=imfilter(I_last,h,'conv','same','replicate');
        
        I_delta=I-I_mean;
        I_delta_mean=median(I_delta(:));
        I_delta_std=mad(I_delta(:),1)*1.4826;
        
        bw_sum_old=sum(bw(:));
        bw=bw|(I_delta>I_delta_mean+I_delta_std*threshold);
        
        bw_sum_diff=sum(bw(:))-bw_sum_old;
        bw=imclose(bw,strel('disk',1,0));
        
        I_last=I;
        I_last(bw)=I_mean(bw)+I_delta_std*threshold;
    end
    
    % Remove single discrete pixels.
    if sum(bw(:))>0
        bwLabel=bwconncomp(bw,ones(3,3));
        area = cellfun(@numel, bwLabel.PixelIdxList);
        idxToKeep=area>=2;
        bwLabel.PixelIdxList=bwLabel.PixelIdxList(idxToKeep);
        bwLabel.NumObjects=sum(idxToKeep);
        
        bw=false(size(I));
        for k = 1 : bwLabel.NumObjects
            bw(bwLabel.PixelIdxList{k}) = true;
        end
    end
end