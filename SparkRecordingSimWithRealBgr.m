function varargout=SparkRecordingSimWithRealBgr(varargin)
    % Simulate confocal recordings of cardiac calcium sparks.
    %
    % This function use the key/value inputs.
    %  
    % [Noisy,Truth,SparkPosition,Bgr]=SparkRecordingSimWithRealBgr(varargin);
    %
    % Input:
    %     KEY           : VALUE   % COMMENT
    %     Gain          :  30
    %     SpatialRes    :  0.28   % In micrometer.
    %     TemporalRes   :  10     % In millisceconds.
    %                   :
    %     Tdim          :  300    % In frames.
    %                   :
    %     SparkNum      :  100
    %                   :
    %     ExpAmp        :  false  % Distribution.
    %     RampAmp       :  false  % Distribution.
    %     SparkAmp      :  1      % In dF/F0.
    %                   :
    %     GaussFWHM     :  false  % Distribution.
    %     CenterFWHM    :  2.8    % In um.
    %     StdFWHM       :  0.5    % In um.
    %                   :
    %     GaussTau      :  false  % Distribution.
    %     CenterTau     :  15     % In ms.
    %     StdTau        :  5
    %                   :
    %     LaserIntensity:  20     % Percentage.
    %     SetupNoise    :  3.5    % a.u.
    %     SetupOffset   :  50     % a.u.
    %                   :
    %     MinTDiff      :  100    % ms.
    %     MinXYDiff     :  5      % um.
    %
    % Output:
    %    SparkWithNoise:    Spark recordings with noise.
    %    SparkNoiseFree:    Spark recordings without noise.
    %    SparkPosition:     (1)ID (2)xc (3)yc (4)t_onset (5)Bgr (6)dF/F0 (7)FWHM/2 (8)Decay.
    %    [SparkMovie,AllSparks,SparkPosition,Bgr]
    
    %% Input.
    p=inputParser;
    p.addParameter('Gain',30,@(x)isscalar(x));
    p.addParameter('SpatialRes',0.28,@(x)isscalar(x));                % In micrometer.
    p.addParameter('TemporalRes',10,@(x)isscalar(x));                 % In millisceconds.
    
    p.addParameter('Tdim',300,@(x)isscalar(x));                       % In frames.
    
    p.addParameter('SparkNum',100,@(x)isscalar(x));
    
    p.addParameter('ExpAmp', false,@(x)islogical(x));                  % Distribution.
    p.addParameter('RampAmp',false,@(x)islogical(x));                  % Distribution.
    p.addParameter('SparkAmp',1,@(x)isscalar(x));                      % In dF/F0.
    p.addParameter('SparkAmpInput',[],@(x)isempty(x) || isvector(x));  % In dF/F0.
    
    p.addParameter('GaussFWHM',false,@(x)islogical(x));                % Distribution.
    p.addParameter('CenterFWHM',2.8,@(x)isscalar(x));                  % In um.
    p.addParameter('StdFWHM',0.5,@(x)isscalar(x));
    p.addParameter('FWHMInput',[],@(x)isempty(x) || isvector(x));      %
    
    p.addParameter('GaussTau', false,@(x)islogical(x));                % Distribution.
    p.addParameter('CenterTau',15,@(x)isscalar(x));                    % In ms.
    p.addParameter('StdTau',   5,@(x)isscalar(x));
    p.addParameter('ExpTau', false,@(x)islogical(x));                  % Distribution.
    p.addParameter('TauInput',[],@(x)isempty(x) || isvector(x));       %
    
    p.addParameter('LaserIntensity',20,@(x)isscalar(x));               % Percentage.
    p.addParameter('SetupNoise',3.5,@(x)isscalar(x));                  % a.u.
    p.addParameter('SetupOffset',50,@(x)isscalar(x));                  % a.u.
    
    p.addParameter('MinTDiff',100,@(x)isscalar(x));                    % ms.
    p.addParameter('MinXYDiff',5,@(x)isscalar(x));                     % um.
    
    p.addParameter('PairDistanceXY',[],@(x)isscalar(x));               % um.
    p.addParameter('PairDistanceT',[],@(x)isscalar(x));                % ms.
    
    p.addParameter('GenerateNoise',true,@(x)islogical(x));             % um.
    
    parse(p, varargin{:});
    p=p.Results;
    
    clear('varargin')
    
    
    Gain=p.Gain;
    SparkNum=p.SparkNum;
    SparkAmpInput=p.SparkAmpInput;
    FWHMInput=p.FWHMInput;
    TauInput=p.TauInput;
    SpatialRes=p.SpatialRes;
    TemporalRes=p.TemporalRes;
    ExpAmp=p.ExpAmp;
    RampAmp=p.RampAmp;
    SparkAmp=p.SparkAmp;
    GaussFWHM=p.GaussFWHM;
    CenterFWHM=p.CenterFWHM;
    StdFWHM=p.StdFWHM;
    GaussTau=p.GaussTau;
    CenterTau=p.CenterTau;
    StdTau=p.StdTau;
    ExpTau=p.ExpTau;
    pTdim=p.Tdim;
    LaserIntensity=p.LaserIntensity;
    SetupOffset=p.SetupOffset;
    SetupNoise=p.SetupNoise;
    MinTDiff=p.MinTDiff;
    MinXYDiff=p.MinXYDiff;
    GenerateNoise=p.GenerateNoise;
    
    PairDistanceXY=round(p.PairDistanceXY/SpatialRes);
    pairedxy=false;
    
    if ~isempty(PairDistanceXY)
        SparkNum=round(SparkNum/2);
        pairedxy=true;
    end
    
    PairDistanceT=round(p.PairDistanceT/TemporalRes);
    pairedt=false;
    if ~isempty(PairDistanceT)
        SparkNum=round(SparkNum/2);
        pairedt=true;
    end
%%
    %% Check gain value.
    if (Gain<=0) || (SparkNum<0)
        disp('    Gain value must be >= 0.')
        disp('    SparkNum value must be >= 0.')
        varargout(1)={[]};
        varargout(2)={[]};
        varargout(3)={[]};
        varargout(4)={[]};
        return
    end
    
    if ~isempty(SparkAmpInput)
        SparkAmpInput=SparkAmpInput(:);
        SparkNum=numel(SparkAmpInput);
    end
    if ~isempty(FWHMInput)
        FWHMInput=FWHMInput(:);
        if SparkNum>numel(FWHMInput)
            SparkNum=numel(FWHMInput);
        end
    end
    if ~isempty(TauInput)
        TauInput=TauInput(:);
        if SparkNum>numel(TauInput)
            SparkNum=numel(TauInput);
        end
    end
    %% Load sample bgr.
    S=SampleCell;
    Bgr=imresize(S.Bgr,S.PixelSize/SpatialRes,'bicubic');
    BgrBW=S.BgrBW;
    BgrBW=imresize(BgrBW,S.PixelSize/SpatialRes,'bicubic');
    if pairedxy
        BgrBWTemp=circshift(BgrBW,[0,-PairDistanceXY]);
        BgrBW=(BgrBW+BgrBWTemp)==2;
        clear('BgrBWTemp')
    end
    
    Bgr=Bgr.*(Bgr>0); %.*imfill(BgrBW,'holes');
    Bgr=Bgr/S.LaserIntensity*LaserIntensity;
    clear('S');
    
    %% Pre-defined parameters.
    RecrdSize=size(Bgr);
    
    %% Pararmeters.
    fprintf('\t==================================================================\n')
    fprintf('\tSpark movie settings:\n')
    fprintf('\t%-30s%0.2f um X %0.2f ms, %0.1f ms/frame\n','Resolution (x,y,t):',SpatialRes,...
        SpatialRes,TemporalRes);
    fprintf('\t%-30s%3.0f x %3.0f x %3.0f\n','Dimension  (x,y,t):',RecrdSize(1),RecrdSize(2),pTdim);
    fprintf('\t%-30s%0.0f%%\n','Laser Intensity:',LaserIntensity);
    fprintf('\t%-30s%0.2f\n','Device Gain:',Gain);
    fprintf('\t%-30s%0.0f\n','Spark number:',SparkNum);
    fprintf('\t%-30s%0.3f\n','Spark amplitude (dF/F0):',SparkAmp);
    fprintf('\t%-30s%0.2f / %0.2f a.u.\n','Detector offset/noise:',SetupOffset,SetupNoise);
    
    
    %% Spark coordinates.
    
    % Shrink the Bgr to avoid edge sparks.
    BgrBW=imerode(BgrBW,strel('disk',round(2/SpatialRes),0));
    
    % Convert Bgr into index.
    BgrIndex=find(BgrBW==1);
    BgrIndexNum=numel(BgrIndex);
    
    % Generate spark position.
    SparkPosIdx=BgrIndex(round(rand(SparkNum,1)*(BgrIndexNum-1))+1);
    [Xcoord,Ycoord]=ind2sub(size(Bgr),SparkPosIdx);
    if pairedt
        Tcoord=floor(rand(SparkNum,1)*(pTdim-1-60/TemporalRes-PairDistanceT))+1;
    else
        Tcoord=floor(rand(SparkNum,1)*(pTdim-1-60/TemporalRes))+1;
    end
    xyt=cat(2,Xcoord,Ycoord,Tcoord);
    clear('Xcoord','Ycoord','Tcoord');
    xyt=sortrows(xyt,3);
    
    for j=1:(SparkNum-1)
        distance=sqrt((xyt(j+1,1)-xyt(j,1))^2+(xyt(j+1,2)-xyt(j,2))^2);
        if ((xyt(j+1,3)-xyt(j,3))<MinTDiff/TemporalRes) && (distance<(MinXYDiff/SpatialRes))
            while distance<(MinXYDiff/SpatialRes)
                SparkPosIdx=BgrIndex(round(rand(1)*BgrIndexNum));
                [Xcoord,Ycoord]=ind2sub(size(Bgr),SparkPosIdx);
                distance=sqrt((Xcoord-xyt(j,1))^2+(Ycoord-xyt(j,2))^2);
            end
            xyt(j+1,1)=Xcoord;
            xyt(j+1,2)=Ycoord;
        end
    end
    if pairedxy
        xyttemp=xyt;
        xyttemp(:,2)=xyttemp(:,2)+PairDistanceXY;
        xyt=cat(1,xyt,xyttemp);
        SparkNum=SparkNum*2;
    end
    if pairedt
        xyttemp=xyt;
        xyttemp(:,3)=xyttemp(:,3)+PairDistanceT;
        xyt=cat(1,xyt,xyttemp);
        SparkNum=SparkNum*2;
    end
    clear('Xcoord','Ycoord','MaskEdge','distance','j','BgrIndex','BgrIndexNum',...
        'SparkPosIdx','xyttemp','pairedxy','pairedt','PairDistanceT','PairDistanceXY');
    
    %% Put spark in.
    AllSparks=zeros(RecrdSize(1),RecrdSize(2),pTdim);
    
    if isempty(SparkAmpInput)
        if ExpAmp
            SaprkAmpDistribution=exprnd(SparkAmp,[SparkNum,1]);
            fprintf('\t%-30sExponential, center = %0.2f\n','Amplitude distribution:',SparkAmp);
        elseif RampAmp
            SaprkAmpDistribution=rand(SparkNum,1)*SparkAmp;
            fprintf('\t%-30sRamp from 0 to %0.2f\n','Amplitude distribution:',SparkAmp);
        else
            fprintf('\t%-30sConstant, center = %0.2f\n','Amplitude distribution:',SparkAmp);
            SaprkAmpDistribution=ones(SparkNum,1)*SparkAmp;
        end
    else
        fprintf('\t%-30sfrom user input\n','Amplitude distribution:');
        SaprkAmpDistribution=SparkAmpInput;
    end
    
    if isempty(FWHMInput)
        if GaussFWHM
            SaprkFWHMDistribution=randn(SparkNum,1)*StdFWHM+CenterFWHM;
            SaprkFWHMDistribution=SaprkFWHMDistribution.*(SaprkFWHMDistribution>0);
            fprintf('\t%-30sGaussian, center/std = %0.2f/%0.2f\n','FWHM distribution:',CenterFWHM,StdFWHM);
        else
            fprintf('\t%-30sConstant, center = %0.2f\n','FWHM distribution:',CenterFWHM);
            SaprkFWHMDistribution=ones(SparkNum,1)*CenterFWHM;
        end
    else
        SaprkFWHMDistribution=FWHMInput;
    end
    SaprkFWHMDistribution=SaprkFWHMDistribution/2;
    
    if isempty(TauInput)
        if GaussTau
            SaprkTauDistribution=randn(SparkNum,1)*StdTau+CenterTau;
            SaprkTauDistribution=SaprkTauDistribution.*(SaprkTauDistribution>0);
            fprintf('\t%-30sGaussian, center/std = %0.2f/%0.2f\n','Tau distribution:',CenterTau,StdTau);
        elseif ExpTau
            SaprkTauDistribution=exprnd(CenterTau,[SparkNum,1]);
            fprintf('\t%-30sExponential, center = %0.2f\n','Tau distribution:',CenterTau);
        else
            SaprkTauDistribution=ones(SparkNum,1)*CenterTau;
            fprintf('\t%-30sConstant, center = %0.2f\n','Tau distribution:',CenterTau);
        end
    else
        SaprkTauDistribution=TauInput;
    end
    fprintf('\t==================================================================\n')
    %%
    SparkPosition=zeros(SparkNum,8);
    for j=1:SparkNum
        CurSpark=SparkDsptSimAllInOneV1_2('Mu',TemporalRes,'SpR',SpatialRes,...
            'TpR',TemporalRes,'OutputSettings',false,...
            'FWHMMax',SaprkFWHMDistribution(j),'Tau',SaprkTauDistribution(j));
        CurrSarkSize=[size(CurSpark,1),size(CurSpark,2),size(CurSpark,3)];
        
        x1=round(xyt(j,1)-CurrSarkSize(1)/2);    x2=x1+CurrSarkSize(1)-1;
        if x1<1;xydiff=1-x1;x1=x1+xydiff;x2=x2+xydiff;end
        if x2>RecrdSize(1);xydiff=x2-RecrdSize(1);x1=x1-xydiff;x2=x2-xydiff;end
        
        y1=round(xyt(j,2)-CurrSarkSize(2)/2);    y2=y1+CurrSarkSize(2)-1;
        if y1<1;xydiff=1-y1;y1=y1+xydiff;y2=y2+xydiff;end
        if y2>RecrdSize(2);xydiff=y2-RecrdSize(2);y1=y1-xydiff;y2=y2-xydiff;end
        
        t1=xyt(j,3);                            t2=t1+CurrSarkSize(3)-1;
        
        CurrBgr=Bgr(x1:x2,y1:y2);
        CurrBgr=median(CurrBgr(:));
        if t2<=pTdim
            AllSparks(x1:x2,y1:y2,t1:t2)=AllSparks(x1:x2,y1:y2,t1:t2)+...
                CurrBgr*CurSpark*SaprkAmpDistribution(j);
        else
            AllSparks(x1:x2,y1:y2,t1:pTdim)=AllSparks(x1:x2,y1:y2,t1:pTdim)+...
                CurrBgr*CurSpark(:,:,1:(pTdim-t1+1))*SaprkAmpDistribution(j);
        end
        SparkPosition(j,:)=[j,xyt(j,2),xyt(j,1),xyt(j,3),CurrBgr,SaprkAmpDistribution(j),...
            SaprkFWHMDistribution(j),SaprkTauDistribution(j)];
    end
    clear('j','t1','t2','x1','x2','y1','y2','xyt','CurrBgr','CurrSarkSize','CurSpark','xydiff')
    
    %% Generate spark movie with cell background.
    % AllSparks=AllSparks.*repmat(BgrBW,[1,1,pTdim]);
    AllSparks=AllSparks+repmat(Bgr,[1,1,pTdim]);
    
    %% Add Poissonian noise.
    if GenerateNoise
        SparkOnCellWithPoissonNoise=zeros(size(AllSparks));
        for k=1:size(AllSparks,3)
            SparkOnCellWithPoissonNoise(:,:,k)=Gain*poissrnd(AllSparks(:,:,k)/Gain);
        end
        % Final movie.
        SparkMovie=SparkOnCellWithPoissonNoise + SetupOffset+...
            randn(RecrdSize(1),RecrdSize(2),pTdim)*SetupNoise;
    else
        SparkMovie=AllSparks;
    end
    %% Add Setup offset without noise.
    AllSparks=AllSparks+SetupOffset;
    
    %% Output
    if nargout==1
        varargout(1)={SparkMovie};
    elseif nargout==2
        varargout(1)={SparkMovie};
        varargout(2)={AllSparks};
    elseif nargout==3
        varargout(1)={SparkMovie};
        varargout(2)={AllSparks};
        varargout(3)={SparkPosition};
        fprintf('\t(1)ID (2)dF/F0 (3)Bgr (4)xc (5)yc (6)t_onset (7)FWHM/2 (8)Decay.\n')
    elseif nargout==4
        varargout(1)={SparkMovie};
        varargout(2)={AllSparks};
        varargout(3)={SparkPosition};
        varargout(4)={Bgr};
        fprintf('\t(1)ID (2)dF/F0 (3)Bgr (4)xc (5)yc (6)t_onset (7)FWHM/2 (8)Decay.\n')
    end
end


