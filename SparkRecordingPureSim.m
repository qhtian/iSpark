function varargout=SparkRecordingPureSim(varargin)
    % Simulate confocal recordings of cardiac calcium sparks.
    % Output:
    %    SparkWithNoise:    Spark recordings with noise.
    %    SparkNoiseFree:    Spark recordings without noise.
    %    SparkPosition:     (1)ID (2)xc (3)yc (4)t_onset (5)Bgr (6)dF/F0 (7)FWHM/2 (8)Decay.
    
    %% Input.
    p=inputParser;
    p.addParameter('Gain',30,@(x)isscalar(x));
    p.addParameter('BasalIntDen',200,@(x)isscalar(x));                 % a.u.
    
    p.addParameter('SetupNoise',3.5,@(x)isscalar(x));                 % a.u.
    p.addParameter('SetupOffset',50,@(x)isscalar(x));                 % a.u.
    p.addParameter('SpatialRes',0.28,@(x)isscalar(x));                % In micrometer.
    p.addParameter('TemporalRes',10,@(x)isscalar(x));                 % In millisceconds.

    p.addParameter('CellArea',2472,@(x)isscalar(x));    
    p.addParameter('Tdim',500,@(x)isscalar(x));                       % In frames.
    
    p.addParameter('SparkNum',100,@(x)isscalar(x));
    p.addParameter('SparkAmp',1,@(x)isscalar(x));                      % In dF/F0.
    p.addParameter('CenterFWHM',2.8,@(x)isscalar(x));                 % In um.
    p.addParameter('StdFWHM',0.5,@(x)isscalar(x));
    p.addParameter('CenterTau',35,@(x)isscalar(x));                   % In ms.
    p.addParameter('StdTau',   5,@(x)isscalar(x));
    
    p.addParameter('MinTDiff',100,@(x)isscalar(x));                   % ms.
    p.addParameter('MinXYDiff',5,@(x)isscalar(x));                    % um.
    
    parse(p, varargin{:});
    p=p.Results;
    clear('varargin')
    
    %% Load sample bgr.
    CellArea=p.CellArea;
    RecordLength=p.Tdim;
    BasalIntDen=p.BasalIntDen;
    Gain=p.Gain;
    GaussSigma=p.SetupNoise;
    Offset=p.SetupOffset;
    SpatialRes=p.SpatialRes;
    TemporalRes=p.TemporalRes;
    SparkNum=p.SparkNum;
    SparkAmp=p.SparkAmp;
    MinTDiff=p.MinTDiff;
    MinXYDiff=p.MinXYDiff;
    CenterFWHM=p.CenterFWHM;
    CenterTau=p.CenterTau;
    
    xLen=sqrt(CellArea/4);
    yLen=xLen*4;
    
    xLen=round(xLen/SpatialRes);
    yLen=round(yLen/SpatialRes);
    
    %% Pararmeters.
    fprintf('\t==================================================================\n')
    fprintf('\tSpark movie settings:\n')
    fprintf('\t%-30s%0.2f um X %0.2f ms X %0.1f ms/frame\n','Resolution (x,y,t):',SpatialRes,SpatialRes,TemporalRes);
    fprintf('\t%-30s%3.0f x %3.0f x %3.0f\n','Dimension  (x,y,t):',xLen,yLen,RecordLength);
    fprintf('\t%-30s%0.2f\n','Device Gain:',Gain);
    fprintf('\t%-30s%0.0f\n','Spark number:',SparkNum);
    fprintf('\t%-30s%0.3f\n','Spark amplitude (dF/F0):',SparkAmp);
    fprintf('\t%-30s%0.2f / %0.2f a.u.\n','Detector offset/noise:',Offset,GaussSigma);
    
    %% Pre-defined parameters.
    BgrBW=false(xLen,yLen+xLen);
    BgrBW(:,1:yLen,:)=true;
    RecrdSize=[xLen,yLen];

    %% Spark coordinates.
    % Convert Bgr into index.
    BgrIndex=find(BgrBW==1);
    BgrIndexNum=numel(BgrIndex);
    
    % Generate spark position.
    SparkPosIdx=BgrIndex(round(rand(SparkNum,1)*(BgrIndexNum-1))+1);
    [Xcoord,Ycoord]=ind2sub(RecrdSize,SparkPosIdx);
    Tcoord=floor(rand(SparkNum,1)*(RecordLength-1-40/TemporalRes))+1;
    xyt=cat(2,Xcoord,Ycoord,Tcoord);
    clear('Xcoord','Ycoord','Tcoord');
    xyt=sortrows(xyt,3);

    for j=1:(SparkNum-1)
        distance=sqrt((xyt(j+1,1)-xyt(j,1))^2+(xyt(j+1,2)-xyt(j,2))^2);
        if ((xyt(j+1,3)-xyt(j,3))<MinTDiff/TemporalRes) && distance<(MinXYDiff/SpatialRes)
            while distance<(MinXYDiff/SpatialRes)
                SparkPosIdx=BgrIndex(round(rand(1)*BgrIndexNum));
                [Xcoord,Ycoord]=ind2sub(RecrdSize,SparkPosIdx);
                distance=sqrt((Xcoord-xyt(j,1))^2+(Ycoord-xyt(j,2))^2);
            end
            xyt(j+1,1)=Xcoord;
            xyt(j+1,2)=Ycoord;
        end
    end
    clear('Xcoord','Ycoord','MaskEdge','distance','j','BgrIndex','BgrIndexNum',...
        'SparkPosIdx');
    
    %% Put spark in.
    AllSparks=zeros(RecrdSize(1),RecrdSize(2)+RecrdSize(1),RecordLength,'single');
    AllSparks(:,1:RecrdSize(2),:)=AllSparks(:,1:RecrdSize(2),:)+BasalIntDen;
    
    fprintf('\t%-30sConstant, center = %0.2f\n','Amplitude distribution:',SparkAmp);
    SaprkAmpDistribution=ones(SparkNum,1)*SparkAmp;
    SaprkFWHMDistribution=ones(SparkNum,1)*CenterFWHM;
    SaprkFWHMDistribution=SaprkFWHMDistribution/2;
    SaprkTauDistribution=ones(SparkNum,1)*CenterTau;
    fprintf('\t%-30sConstant, center = %0.2f\n','Tau distribution:',CenterTau);

    fprintf('\t==================================================================\n')
    %%
    SparkPosition=zeros(SparkNum,8);
    for j=1:SparkNum
        CurSpark=SparkDsptSimAllInOneV1_2('Mu',TemporalRes,'SpR',SpatialRes,...
            'TpR',TemporalRes,'OutputSettings',false,...
            'FWHMMax',SaprkFWHMDistribution(j),'Tau',SaprkTauDistribution(j));
        CurrSarkSize=[size(CurSpark,1),size(CurSpark,2),size(CurSpark,3)];
        
        x1=round(xyt(j,1)-CurrSarkSize(1)/2);
        x2=x1+CurrSarkSize(1)-1;
        if x1<1
            xydiff=1-x1;
            x1=x1+xydiff;
            x2=x2+xydiff;
        end
        if x2>RecrdSize(1)
            xydiff=x2-RecrdSize(1);
            x1=x1-xydiff;
            x2=x2-xydiff;
        end
        
        y1=round(xyt(j,2)-CurrSarkSize(2)/2);
        y2=y1+CurrSarkSize(2)-1;
        if y1<1
            xydiff=1-y1;
            y1=y1+xydiff;
            y2=y2+xydiff;
        end
        if y2>RecrdSize(2)
            xydiff=y2-RecrdSize(2);
            y1=y1-xydiff;
            y2=y2-xydiff;
        end
        
        t1=xyt(j,3);
        t2=t1+CurrSarkSize(3)-1;
        
        CurrBgr=BasalIntDen;
        if t2<=RecordLength
            AllSparks(x1:x2,y1:y2,t1:t2)=AllSparks(x1:x2,y1:y2,t1:t2)+...
                CurrBgr*CurSpark*SaprkAmpDistribution(j);
        else
            AllSparks(x1:x2,y1:y2,t1:RecordLength)=AllSparks(x1:x2,y1:y2,t1:RecordLength)+...
                CurrBgr*CurSpark(:,:,1:(RecordLength-t1+1))*SaprkAmpDistribution(j);
        end
        SparkPosition(j,:)=[j,xyt(j,2),xyt(j,1),xyt(j,3),CurrBgr,SaprkAmpDistribution(j),SaprkFWHMDistribution(j),SaprkTauDistribution(j)];
    end
    clear('j','t1','t2','x1','x2','y1','y2','xyt','CurrBgr','CurrSarkSize','CurSpark','xydiff')
    
    %% Generate spark movie with cell background.
    AllSparks=bsxfun(@times,AllSparks,BgrBW);
    
    %% Add Poissonian noise.
    SparkOnCellWithPoissonNoise=zeros(size(AllSparks),'like',AllSparks);
    for k=1:size(AllSparks,3)
        SparkOnCellWithPoissonNoise(:,:,k)=Gain*poissrnd(AllSparks(:,:,k)/Gain);
    end
    
    % Final movie.
    AllSparks=AllSparks+Offset;
    SparkMovie=SparkOnCellWithPoissonNoise+Offset+randn(size(SparkOnCellWithPoissonNoise))*GaussSigma;
    
    %% Output
    varargout(1)={SparkMovie};
    varargout(2)={AllSparks};
    varargout(3)={SparkPosition};
    fprintf('\t(1)ID (2)dF/F0 (3)Bgr (4)xc (5)yc (6)t_onset (7)FWHM/2 (8)Decay.\n')

end


