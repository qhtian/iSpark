function [SparkPos,SparkProperty,allFitting]=SingleSparkv1_2_Parallel(raw,...
          sparkprofile,SparkLabel,SparkTLim,xyt_dim,CellMask)
%% SingleSparkv1_2: check and extract information for every single spark.
%
%  ======================================
%  Just for internal usage. Qinghai Tian.
%  ====================================== 
%
%  SparkPos=[x1,x2,y1,y2,z1,z2,xc,yc,z_onset,No.];
%  SparkProperty=[FWHM,MajorAxisLength,MinorAxisLength,Amp,Bgr,Amp/Bgr,Tau,
%                 Sigma,DetctionMass,dF/F0_Mass,SpToMembraneDistance,No.];

%% Presettings.
    %pixel, extend the spatial dimension to get reliable FWHM and better amplitude.
    SparkMaskExtension=ceil(1.68/xyt_dim(1));

    %[leading foot length, tail length]
    SparkTLim=round(SparkTLim/xyt_dim(3));
    spark_num=numel(sparkprofile);
    RawSize1=size(raw,1);RawSize2=size(raw,2);RawSize3=size(raw,3);
    
    SparkPosExtended=zeros(spark_num,9);

    SparkPos=zeros(spark_num,9);
    SpCenterX=zeros(spark_num,1);
    SpCenterY=zeros(spark_num,1);
    SparkProperty=zeros(spark_num,12);

    rawsmthn_RawData=cell(spark_num,1);
    sparkmask_RawData=cell(spark_num,1);
    sparkmaskOther_RawData=cell(spark_num,1);
    % RoughBgr_RawData=cell(spark_num,1);
    frameFittedData=cell(spark_num,1);
    gauWxpFittedData=cell(spark_num,1);
    
    GauExpFit_InitialMu=zeros(spark_num,1);
    GauExpFit_RawData=cell(spark_num,1);
    % BgrFrame_Bgr=zeros(spark_num,1);
    SpMss=zeros(spark_num,1);
    
    z11_All=zeros(spark_num,1);
    sparkprofile_Area=zeros(spark_num,1);
    % sparkprofile_Orientation=zeros(spark_num,1);
    
    RoughBgr=median(raw,3);
    raw_smthn=raw-repmat(RoughBgr,[1,1,size(raw,3)]);

%% Cycle single spark.
    for k=1:spark_num
        % Get the detected region.
        z1=ceil(sparkprofile(k,1).BoundingBox(3));
        z2=floor(sparkprofile(k,1).BoundingBox(3)+sparkprofile(k,1).BoundingBox(6)); 
        x1=ceil(sparkprofile(k,1).BoundingBox(2));
        x2=floor(sparkprofile(k,1).BoundingBox(2)+sparkprofile(k,1).BoundingBox(5));         
        y1=ceil(sparkprofile(k,1).BoundingBox(1));
        y2=floor(sparkprofile(k,1).BoundingBox(1)+sparkprofile(k,1).BoundingBox(4)); 
        SparkPos(k,1:6)=[x1,x2,y1,y2,z1,z2];

        % Extend the detected region.
        x1=x1-SparkMaskExtension;        if x1<1;x1=1;end;if x1>RawSize1;x1=RawSize1;end
        x2=x2+SparkMaskExtension;        if x2<1;x2=1;end;if x2>RawSize1;x2=RawSize1;end
        y1=y1-SparkMaskExtension;        if y1<1;y1=1;end;if y1>RawSize2;y1=RawSize2;end
        y2=y2+SparkMaskExtension;        if y2<1;y2=1;end;if y2>RawSize2;y2=RawSize2;end
        if z1<1;z1=1;end;if z1>RawSize3;z1=RawSize3;end
        if z2<1;z2=1;end;if z2>RawSize3;z2=RawSize3;end
        SparkPosExtended(k,1:6)=[x1,x2,y1,y2,z1,z2];
 
        
        % Find out the spatial center of the spark.
        rawsmthn_curr=sum(raw_smthn(x1:x2,y1:y2,z1:z2),3); %#ok<*PFBNS>
        
        MassRawMask=SparkLabel(x1:x2,y1:y2,z1:z2)==k;
        sparkmask_curr=sum(MassRawMask,3)>0;
        SpMss(k)=sum(MassRawMask(:))*xyt_dim(3)*xyt_dim(1)*xyt_dim(2);
        clear('MassRawMask')

        M=imresize(rawsmthn_curr,10,'bicubic').*imresize(sparkmask_curr,10,'nearest');
        [~,Mxy]=max(M(:));
        [cx_curr,cy_curr]=ind2sub(size(rawsmthn_curr)*10,Mxy);
        cx_curr=mean(cx_curr)/10;    cy_curr=mean(cy_curr)/10;

        SparkPosExtended(k,7:8)=[cx_curr,cy_curr];
        clear('M','rawsmthn_curr','Mxy')

        % Save the center of the spark.
        SpCenterX(k)=cx_curr+x1-1;
        SpCenterY(k)=cy_curr+y1-1;
        SparkPos(k,7)=SpCenterX(k);
        SparkPos(k,8)=SpCenterY(k);

        % Find temporal center
        x1_for_t=floor(SpCenterX(k));	if x1_for_t<1;x1_for_t=1;end
        x2_for_t=ceil(SpCenterX(k));	if x2_for_t>RawSize1;x2_for_t=RawSize1;end
        y1_for_t=floor(SpCenterY(k));	if y1_for_t<1;y1_for_t=1;end
        y2_for_t=ceil(SpCenterY(k));	if y2_for_t>RawSize2;y2_for_t=RawSize2;end
        M=raw_smthn(x1_for_t:x2_for_t,y1_for_t:y2_for_t,z1:z2);
        M=sum(sum(M,1),2);  M=M(:);
        cz_curr=find(M==max(M));cz_curr=cz_curr(1);
        SpCenterZ=cz_curr+z1-1;
        
        SparkPos(k,9)=SpCenterZ;
        SparkPosExtended(k,9)=SpCenterZ;
        clear('M','x1_for_t','x2_for_t','y1_for_t','y2_for_t','cz_curr','cx_curr','cy_curr')
        

        % FWHM calculation with frame fitting.
        sparkmask_other=mean((SparkLabel(x1:x2,y1:y2,z1:z2)~=k) & ...
                            (SparkLabel(x1:x2,y1:y2,z1:z2)>0),3)>0;
        rawsmthn_curr=mean(raw_smthn(x1:x2,y1:y2,z1:z2),3);
        
        rawsmthn_RawData{k,1}=rawsmthn_curr;
        sparkmask_RawData{k,1}=sparkmask_curr;
        sparkmaskOther_RawData{k,1}=sparkmask_other;
        % RoughBgr_RawData{k,1}=RoughBgr(x1:x2,y1:y2);
        clear('sparkmask_curr','rawsmthn_curr','sparkmask_other')

        % Spark transient fitting.
        % Extend in temporal dimension.
        Mu=z1;      z1=z1-SparkTLim(1); if z1<1;z1=1;end;   if z1>RawSize3;z1=RawSize3;end
        Mu=Mu-z1+1; z2=z2+SparkTLim(2); if z2<1;z2=1;end;   if z2>RawSize3;z2=RawSize3;end
        x1=floor(SpCenterX(k)); x2=ceil(SpCenterX(k));
        y1=floor(SpCenterY(k)); y2=ceil(SpCenterY(k));
        if x1<1;x1=1;end;if x1>RawSize1;x1=RawSize1;end
        if x2<1;x2=1;end;if x2>RawSize1;x2=RawSize1;end
        if y1<1;y1=1;end;if y1>RawSize2;y1=RawSize2;end
        if y2<1;y2=1;end;if y2>RawSize2;y2=RawSize2;end
        
        SparkPosExtended(k,5)=z1;

        M=raw(x1:x2,y1:y2,z1:z2);

        % Get a safe range of spark transient.
        Mask=(SparkLabel(x1:x2,y1:y2,z1:z2)~=k) & (SparkLabel(x1:x2,y1:y2,z1:z2)>0);
        Mask=sum(sum(Mask,1),2);    Mask=Mask(:)>0;
        cz_curr=(floor(SpCenterZ)-z1+1);
        m=cz_curr;
        while m>=1 && ~Mask(m)
            m=m-1;
        end
        z11=m+1;
        z11_All(k)=z11;
        
        Mu=Mu-z11+1;
        GauExpFit_InitialMu(k)=Mu;
        
        m=cz_curr; while m<=numel(Mask) && ~Mask(m); m=m+1; end; z2=m-1;
        
        M=M(:,:,z11:z2);

        % Get central pixel transient. New Codes.
        M=calCenterTransient(M,SpCenterX(k),SpCenterY(k));
        
        GauExpFit_RawData{k}=M;
        sparkprofile_Area(k)=sparkprofile(k,1).Area*xyt_dim(3)*xyt_dim(1)*xyt_dim(2);
        
        
        clear('Mu','x1','x2','y1','y2','z1','z2','SpCenterZ','Mask','cz_curr','M','m','z11');
    end

    clear('M','RawSize1','RawSize2','RawSize3','SparkLabel','SparkTLim','k','raw','raw_smthn','sparkprofile');


%% parallel fitting.    
    parfor k=1:spark_num
        % Frame fitting.
        [FWHM,disI_pair]=PeakFrameFitting(rawsmthn_RawData{k,1},...
            sparkmask_RawData{k,1},sparkmaskOther_RawData{k,1},...
            SparkPosExtended(k,7),SparkPosExtended(k,8),SparkMaskExtension);
        FWHM=FWHM*xyt_dim(1)*2;
        frameFittedData{k,1}=disI_pair;
        
        sparkMorphology=regionprops(sparkmask_RawData{k,1},'Eccentricity','MajorAxisLength','MinorAxisLength','Orientation');
        Axis_Max=sparkMorphology(1).MajorAxisLength*xyt_dim(1);
        Axis_Min=sparkMorphology(1).MinorAxisLength*xyt_dim(1);

        % Spark mass.
        % BgrFrame=BgrFrame+mean(RoughBgr_RawData{k,1}(:));
        
        % if BgrFrame~=0
        %     SpMss_curr=SpMss(k)/BgrFrame;
        % else
        %     SpMss_curr=NaN;
        % end

        % Transient fitting.
        [Amp,Bgr,MuFit,Tau,Sigma,~,gauExpFit]=GauExpFit(GauExpFit_RawData{k},GauExpFit_InitialMu(k));
        gauWxpFittedData{k,1}=gauExpFit;
        Tau=Tau*xyt_dim(3);
        Sigma=Sigma*xyt_dim(3);
        
        SpCenterZ=MuFit+z11_All(k)-2+SparkPosExtended(k,5);
        SparkPos(k,9)=SpCenterZ;
        
        sp2memdist=SparkToMembraneDistance([ SpCenterX(k),SpCenterY(k)],CellMask)*xyt_dim(1);
        
        SparkProperty(k,:)=[FWHM,Axis_Max,Axis_Min,Amp,Bgr,Amp/Bgr,Tau,Sigma,...
            sparkprofile_Area(k),SpMss(k),sp2memdist,sparkMorphology(1).Orientation];
    end

%% Add spark ID.
    SparkPos=cat(2,SparkPos,(1:spark_num)');
    SparkProperty=cat(2,SparkProperty,(1:spark_num)');
    allFitting=cat(2,frameFittedData,gauWxpFittedData,rawsmthn_RawData,sparkmask_RawData);
end

function dist=SparkToMembraneDistance(Pos,Mask)
    Mask=(Mask-ordfilt2(Mask, 1, [0 1 0;1 1 1;0 1 0],'zeros'))==1;
    [X,Y]=ndgrid(1:size(Mask,1),1:size(Mask,2));
    X=X(Mask);Y=Y(Mask);
    dist=min(sqrt((Pos(1)-X).^2+(Pos(2)-Y).^2));
end

%%
function [FWHM,disI_pair]=PeakFrameFitting(I,mask,maskother,cx,cy,ext)
    I=double(I);
    [sizey,sizex] = size(I);
%% Get center of mass, amplitude, and sigma.
    [X,Y]=ndgrid(1:sizey,1:sizex);
    distance=sqrt((X-cx).^2+(Y-cy).^2); 
    sigma=sqrt(sum(mask(:))/pi)/2;

    mask=imdilate(mask,strel('disk',ext, 0));
    mask=mask & (~maskother);
    distance=distance(mask);    I=I(mask);
    disI_pair=cat(2,distance,I);
    disI_pair=sortrows(disI_pair,1);

    I_max=max(I);
    Dis_max=max(distance);

%% Do a Gaussian fitting with mu=0 and Bgr=0.
    Gau=@(x,p)p(1)*exp(-x.^2/2/p(2)^2)+p(3);
    fun_dev=@(x,y,p)sum((Gau(x,p)-y).^2);
    options=optimset('MaxIter',100000000,'Display','off');
    p0=[I_max sigma 0];
    pLB=[I_max/10 sigma/10 -abs(I_max*3)];
    pUB=[I_max*20 Dis_max*3 I_max];
    if (numel(p0)~=numel(pLB)) || (numel(p0)~=numel(pLB))
        FWHM=NaN;
        disI_pair=NaN(1,2);
    else
        p=fminsearchbnd(@(p)fun_dev(disI_pair(:,1),disI_pair(:,2),p),p0,pLB,pUB,options);
        
        fitted=Gau(disI_pair(:,1),p);
        % Bgr=min(fitted);
        
        disI_pair=cat(2,disI_pair,fitted);
        %% Final results;
        % R2=1-sum((disI_pair(:,2)-fitted).^2)/sum((disI_pair(:,2)-mean(I)).^2);
        FWHM=p(2);
    end
end



%% GauExpFit
function [Amp,Bgr,Mu,Tau,Sigma,R2,fitted]=GauExpFit(Ca,Mu)
% GauConvExp=@(x,AmpMax,Bgr, Mu,Tau,Sigma) AmpMax/2.*exp((Sigma^2+2*Tau*(Mu-...
%              x))/2/Tau^2).*(1-(Sigma^2+Tau*(Mu-x))./abs(Sigma^2+Tau*(Mu-...
%              x)).*erf(abs(Sigma^2+Tau*(Mu-x))/sqrt(2)/Sigma/Tau))+Bgr;

%% Input.
    Ca=double(Ca);
    minCa=min(Ca);
    maxCa=max(Ca);
    numCa=numel(Ca);
    
    if numel(Ca)<3
        Amp=-1; Bgr=-1; Mu=0; Tau=-1; Sigma=-1; R2=0; fitted=[];
        return
    end
    x=(1:numCa)';
    if Mu<1
        Mu=find(Ca==maxCa);  % Mu=Mu(1)-10/tdim;
        if numel(Mu)>1
            Mu=Mu(1);
        end
    end

%% Fit the temporal parameters firstly.
    GauConvExp1=@(x,p) p(4)/2.*exp((p(3)^2+2*p(2)*(p(1)-x))/2/p(2)^2).*(1-(p(3)^2+...
                       p(2)*(p(1)-x))./abs(p(3)^2+p(2)*(p(1)-x)).*erf(abs(p(3)^2+p(2)*(p(1)-...
                       x))/sqrt(2)/p(3)/p(2)))+p(5);
    fun_dev1=@(x,y,p)sum((GauConvExp1(x,p)-y).^2);
    options=optimset('MaxIter',100000000,'Display','off','TolX',1e-4);

    Tau=sum((Ca-minCa)>0.6*(maxCa-minCa));    % if Tau<0; Tau=15/tdim; end
    Sigma=Tau/5;
    
    Bgr=sort(Ca);   Bgr=Bgr(max(1,round(numel(Bgr)*0.15))); 
    Amp=(maxCa-Bgr)/0.6;
    
    p0=[Mu,Tau,Sigma,Amp,Bgr];
    pLB=[0,Tau/10,0,0,minCa];
    pUB=[numCa/2,Tau*3,numCa/2,(maxCa-minCa)*2,median(Ca)];
    
    if (numel(p0)~=numel(pLB)) || (numel(p0)~=numel(pLB))
        Amp=NaN;
        Bgr=NaN;
        Mu=NaN;
        Tau=NaN;
        Sigma=NaN;
        R2=0;
        fitted=NaN;
    else
        p=fminsearchbnd(@(p)fun_dev1(x,Ca,p),p0,pLB,pUB,options);
        Mu=p(1);
        Tau=p(2);
        Sigma=p(3);
        fitted=GauConvExp1(x,p);
        Bgr=min(fitted);
        Amp=max(fitted)-Bgr;
        R2=1-sum((Ca-fitted).^2)/sum((Ca-mean(Ca)).^2);
        
        fitted=cat(2,x,Ca,fitted);
    end
end

function M=calCenterTransient(M,x,y)
    if size(M,1)>1
        x1=floor(x);    x2=ceil(x);
        M=(1-(x-x1))*M(1,:,:)+(1-(x2-x))*M(2,:,:);
    end
    
    if size(M,2)>1
        y1=floor(y);    y2=ceil(y);
        M=(1-(y-y1))*M(:,1,:)+(1-(y2-y))*M(:,2,:);
    end
    M=M(:);
end