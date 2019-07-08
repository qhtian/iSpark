function varargout=SparkAnalysis(RawImage,varargin)
    % An analysis algorithm for the cardiac calcium release events, the calcium sparks.
    %
    % Matlab Image Processing Toolbox, Paralle Computing Toolbox and Statistics
    %  Toolbox are required to run this program. This Program was validated with
    %  Matlab 2012a, and not tested with other versions.
    %
    % Cell Mask Detection: Cell mask detection was based on the isodata method.
    % See paper from T.W. Ridler, S. Calvard, Picture thresholding using an
    % iterative selection method, IEEE Trans. System, Man and Cybernetics,
    % vol.8, 630-632, 1978.
    %
    % Poissonian Noise Handling: Adaptive iterative denoising. See paper from
    % F. Luisier, C. Vonesch, T. Blu, M. Unser, Fast Interscale Wavelet Denoising
    % of Poisson-corrupted Images, Signal Processing, vol. 90, pp. 415-427, 2010.
    % An alternative is ImageJ PURE-denoise plugin developed by the authors.
    % Gaussian noise Handling: discrete consine transform-based GCV minimization.
    % For minizaing GCV score, see paper from Garcia D, Robust smoothing of gridded
    % data in one and higher dimensions with missing values. Computational
    % Statistics & Data Analysis, 2010.
    %
    %
    % Usage of the main program:
    %   Result=SparkAnalysis(ImageStack,varargin);
    %
    %   varargin syntax:
    %       1.  'xyt_dim',[0.28 0.28 10]; Spatial x/y dimension 0.56/0.56 um and
    %            temporal dimension 10 ms/frame.
    %       2.  'Sensitivity', Less false positives (1), sensitive (2),  very
    %            sensitive (3), and very very sensitive (4);
    %       3.  'Threshold',[]; The pixel signal larger than the 'Threshold'
    %            times of noise sigma will be treated as a pixel of sparks.
    %            If it is left empty, the program will determin it automatically.
    %       4.  'DetectionScale', 20~60ms. The rough temporal length to cover typical
    %            sparks.
    %       5.  'GauAmpLim',10; The maximum limit for the fitted amplitude of a spark
    %            in F/F0.
    %       6.  'TauLim', 500, The maximum limit (ms) for the fitted decay constant.
    %       7.  'DetectionLimit',[6.27 60000]; The volume limits of spark detection in
    %            squm*ms. The physical limits are [0.56,0.56,20] ~ [10,10,600] (um,um,ms).
    %       8.  'FWHMLimit',[0.2 10] um; The spark FWHM limits.
    %       9.  'SparkLim',[50 100] ms; A rough range in pixles for a spark transient,
    %            in a format like [leading foot length, tail length];
    %       10. 'PUREDenoising', 0 (skip), 1 (2D) or 2 (3D); The process of PURE-LET denoising.
    %       11. 'GaussDenoise', 0; 0 means GCV minimization with smoothn, positive number
    %		     means manual parameter s for smoothn. A negative number means skip the
    %            denoising process.
    %       12.  'MinSparkSiteDistance', 1. The minimum distance between two
    %            separated spark sites.
    %
    %   Result has the following fields:
    %       1.  DetectionLimit:     the same parameter from the input.
    %       2.  GauAmpLim:   the same parameter from the input.
    %       3.  Threshold:   the same parameter from the input.
    %       4.  SparkLim:    the same parameter from the input.
    %       5.  xyt_dim:     the same parameter from the input.
    %       6.  Creator:     the program used to generate the results.
    %       7.  Data:        the denoised data of input image stack.
    %       8.  Gain:        the original device gain value.
    %       9.  CellMask:    detected cell mask.
    %       10. CellMaskEdge: detected cell edge from CellMask.
    %       11. CameraOffset: dark background of the detector.
    %       12. GaussNoise:   sigma of dark noise.
    %       13. SparkLabel:   all the sparks detected and labeled in specific number.
    %       14. SparkPos: the coordinates of the final sparks, in the format of
    %                         [x1,x2,y1,y2,t1,t2,xc,yc,t_onset,ID_In_SparkLabel].
    %       15. SparkProperty: the calculated spark properties, in the format of
    %                         [FWHM,FWHM_R2,Amp,Bgr,Tau,Sigma,Amp_R2,Old_ID];
    %       16. SparkSitePos: same information as SparkPos, but for spark sites.
    %       17. SparkSiteProperty: same information as SparkProperty, but for spark sites.
    %       18. Spark_Site_Relation: the mapping of single sparks to spark sites.
    %       19. MaximaMask_1_Raw: Detected spark mask in first round.
    %       20. SparkPosRejected: the coordinates of the rejected sparks, in the format
    %                         of [x1,x2,y1,y2,t1,t2,xc,yc,t_onset,ID_In_SparkLabel].
    %       21. SparkPropertyRejected: the rejected spark properties, in the format of
    %                         [FWHM,FWHM_R2,Amp,Bgr,Tau,Sigma,Amp_R2,Old_ID,AmpCheck,
    %                          FWHMCheck,MassCheck];
    %       22. MaximaMask_2_CheckingKept: finally detected spark masks.
    %       23. MaximaMask_2_CheckingRejected: detected but finally rejected spark masks.
    %       24. SparkFrequency: spark frequency and spark site frequency (/squm/s).
    %
    %
    % Notes:
    %   1. Spark transient is fitted with Gau*Exp function:
    %       GauConvExp=@(x,AmpMax,Bgr, Mu,Tau,Sigma)
    %                   AmpMax/2.*exp((Sigma^2+2*Tau*(Mu-x))/2/Tau^2).*(1-...
    %                   (Sigma^2+Tau*(Mu-x))./abs(Sigma^2+Tau*(Mu-x)).*erf(...
    %                   abs(Sigma^2+Tau*(Mu-x))/sqrt(2)/Sigma/Tau))+Bgr;
    %       AmpMax:     maximum amplitude of the spark.
    %       Bgr:        baseline of the spark transient in the center of the spark.
    %       Mu:         the time point of maximum Ca release.
    %       Tau:        Ca removal constant.
    %       Sigma:      the FWHM duration of the Ca release of the spark.
    

    %% Save the programe version.
    Creator='iSpark_V5.02_GUI';
    
    %% Start counting time.
    totaltime=tic;
    
    %% Pruse inputs.
    RawSize=size(RawImage);
    S=inputParser;
    S.addParameter('xyt_dim',       [0.28,0.28,10],@(x)numel(x)==3 && min(x)>0);         % um, um, ms/frame
    S.addParameter('Threshold',[],                 @(x)(isscalar(x) && x>0) || isempty(x));
    S.addParameter('DetectionScale',    [20,10,60],@(x)numel(x)==3 && min(x)>=1);        % in ms.
    S.addParameter('GauAmpLim',                 10,@(x)isscalar(x) && x>0);
    S.addParameter('CameraOffset',             NaN,@(x)isscalar(x) && ((x>=0)|| (isnan(x))));
    S.addParameter('TauLim',                  1500,@(x)isscalar(x) && x>0);
    S.addParameter('DetectionLimit',  [6.27,60000],@(x)numel(x)==2 && min(x)>=0);        % Detection mass limit.
    S.addParameter('FWHMLimit',           [0.2,10],@(x)numel(x)==2 && min(x)>0);         % Min/Max FWHM of spark.
    S.addParameter('SparkLim',            [50,150],@(x)numel(x)==2 && min(x)>0)          % [leading foot length, tail length]
    S.addParameter('Denoising',                  1,@(x)numel(x)==1 && (x==1 || x==2 || x==0));  % 1: PUREdenoise! 0: skip; 2: CANDLEdenoise
    S.addParameter('CellMasking',                1,@(x) ((size(x,1)==RawSize(1)) && ...
        (size(x,2)==RawSize(2))) ||(x==1) || (x==2) || (x==3) || (x==4));
        % 1: Isodata; 2. Local adaptive threshold; 3. Fine Masking; 4. User Input.
    
    S.addParameter('PhaseShiftCorrect',         0,@(x)((numel(x)==1) && ((x==0)||(x==1)||(x==2))));
    S.addParameter('MinSparkSiteDistance',      1,@(x)isscalar(x) && x>0)               % In um.
    S.addParameter('ImRegister',            false,@(x)islogical(x));
    S.addParameter('ApplyGlobalCellmask',    true,@(x)islogical(x));
                % Apply the detected mask on the global spark mask.
                % Important for PMT/HyD-based recordings to remove
                % spike pixels from the background.
    S.addParameter('UseWatershedSegment',    true,@(x)islogical(x));
    S.addParameter('ApplyImdilation',       false,@(x)islogical(x));
    S.addParameter('PhotonCountingFlow',       [],@(x) isempty(x) || ((size(x,1)==RawSize(1)) && (size(x,2)==RawSize(2))));
    parse(S, varargin{:});
    S=S.Results;
    clear('varargin')
    S.Creator=Creator; clear('Creator')

    %% Hard cut off parameters.
    LowestPixelPhotonCount=0.5;
    
    %% Sort parameters to avoid confusings.
    S.DetectionLimit=sort(S.DetectionLimit,'ascend');
    S.FWHMLimit=sort(S.FWHMLimit,'ascend');

    %% Ignore pixel dimension difference bewteen x and y.
    xy_diff=S.xyt_dim(1)/S.xyt_dim(2);
    if xy_diff>1.01 || xy_diff<0.99
        fprintf('\n    ---------------------------- WARNING ----------------------------\n');
        fprintf(['    Pixel dimension difference between X and Y exceeds 1.0%% ' ...
            '(%0.1f%%).\n    iSpark will ignore it. You should take care of it.\n'],...
            abs(xy_diff-1)*100);
        fprintf('    ---------------------------- WARNING ----------------------------\n\n');
    end
    S.xyt_dim=[(S.xyt_dim(1)+S.xyt_dim(2))/2,(S.xyt_dim(1)+S.xyt_dim(2))/2,S.xyt_dim(3)];
    clear('xy_diff');

    %% Display predefined parameters.
    fprintf('    Spark Analysis (Ver. %s)\n',S.Creator);
    disp('    =======================================================================');
    disp('    Predefined settings:')
    
    Str='    Physical dimensions:';
    fprintf('%-50s%0.2fum x %0.2fum x %0.2fms\n',Str,S.xyt_dim(1),S.xyt_dim(2),S.xyt_dim(3));
    
    Str='    Threshold strategy :';
    if ~isempty(S.Threshold)
        fprintf('%-50sthreshold = %0.2f\n',Str,S.Threshold);
    else
        fprintf('%-50sADAPTIVE Threshold\n',Str);
    end
    
    Str='    Detection scale    :';
    fprintf('%-50s[%0.0f, %0.0f] ms\n',Str,S.DetectionScale(1),S.DetectionScale(3));
    
    Str='    Maximum amplitude  :';
    fprintf('%-50sdF/F0 <= %1.1f\n',Str,S.GauAmpLim);
    
    Str='    Min/Max spark mass :';
    fprintf('%-50s[%g, %g] squm*ms\n',Str,S.DetectionLimit(1),S.DetectionLimit(2));
    
    Str='    Min/Max decay Tau  :';
    fprintf('%-50s[%g, %g] ms\n',Str,0,S.TauLim);
    
    Str='    Min/Max FWHM       :';
    fprintf('%-50s[%g, %g] um\n',Str,S.FWHMLimit(1),S.FWHMLimit(2));
    
    Str='    Leading Foot/Tail  :';
    fprintf('%-50s[%g, %g] ms\n',Str,S.SparkLim(1),S.SparkLim(2));
    
    disp('    =======================================================================');
    clear('Str','str1');

    %% Put the raw data into the struct I.
    S.Data=single(RawImage);    clear('RawImage');
    
    Str='    Input data size:';
    fprintf('%-50s%g x %g x %g\n',Str,RawSize(1),RawSize(2),RawSize(3));
    clear('Str');

    %% Phase contrast correction.
    if S.PhaseShiftCorrect>0
        fprintf('%-50s','    Phase contrast correcting:');
        [S.Data,S.PhaseShiftCorrect]=PhaseShiftAutoAlign(S.Data,S.PhaseShiftCorrect);
        fprintf('shift  = %0.3f pixel\n',-S.PhaseShiftCorrect);
    end
    
    %% Cell mask detection
    if numel(S.CellMasking)==1
        [S.CellMask,algo]=CellMasking(S.Data,S.CellMasking,S.xyt_dim(1:2));
    else
        S.CellMask=S.CellMasking;
        S.CellMasking=4;
        algo='manual checked mask';
    end
    Str=['    Cell masked with ',algo,':'];
    fprintf('%-50sarea   = %0.1f squm\n',Str,S.xyt_dim(1)*S.xyt_dim(2)*sum(S.CellMask(:)));
    clear('str','algo')

    S.CellMaskEdge=(S.CellMask-ordfilt2(S.CellMask, 1, [0,1,0;1,1,1;0,1,0],'zeros'))==1;
    
    %% Denoise based on the estimated gain value.
    externalOffset=S.CameraOffset;
    if S.Denoising==2
        fprintf('%-50s','    CANDLE Poissonian/Gaussian denoising:');
        [S.CameraOffset,S.Gain]=deviceOffset(S.Data);
        S.GaussNoise=NaN;
        S.Data=CANDLEdenoise(S.Data);
        fprintf('done\n');
    elseif S.Denoising==1
        fprintf('%-50s','    PURE-LET Poissonian/Gaussian denoising:');
        [S.Data,S.Gain,S.CameraOffset,S.GaussNoise]=PUREDenoiseJava(S.Data,'CS',6);
        fprintf('done\n');
    else
        S.GaussNoise=NaN;
        S.Gain=NaN;
    end
    if ~isnan(externalOffset)
        S.CameraOffset=externalOffset;
    elseif S.Denoising==0
        [S.CameraOffset,S.Gain]=deviceOffset(S.Data);
    end
    S.Data=S.Data-S.CameraOffset;
    clear('externalOffset');

    Str='    Detector Gain calculated:';
    fprintf('%-50sgain   = %0.2f a.u.\n',Str,S.Gain);
    Str='    Detector Offset calculated:';
    fprintf('%-50soffset = %0.2f a.u.\n',Str,S.CameraOffset);
    Str='    Detector Thermal Noise:';
    fprintf('%-50ssigma  = %0.2f a.u.\n',Str,S.GaussNoise);

    %% Image shift correction.
    if S.ImRegister
        fprintf('%-50s','    Image translation correction...');
        S.Data=ImAutoRegister(S.Data);
        fprintf('done\n');
    end
    
    %% Local spark thresholding.
    Str='    Thresholding sparks:';
    fprintf('%-50s',Str);
    ScaleJ=round(S.DetectionScale/S.xyt_dim(3));
    ScaleJ(ScaleJ<1)=1;
    if ScaleJ(1)<2
        ScaleJ(1)=2;
    end

    DataCV=mad(diff(S.Data,1,3),1,3)*1.4826;
    DataCV=medfilt2(DataCV,[3,3],'symmetric');
    DataCV(DataCV<0)=0;
    DataIntensity=mean(S.Data,3);
    DataIntensity(DataIntensity<0)=0;
    DataCV=DataCV./DataIntensity;
    DataCV(isnan(DataCV))=inf;
    S.CVofData=100*median(DataCV(S.CellMask));
    
    if isempty(S.Threshold)
        S.Threshold=AdaptiveThreshold(DataCV,S.Denoising);
    end
    
    [MaximaMask,S.Threshold]=SparkThresholding(S.Data,'ScaleN',ScaleJ,'Threshold',S.Threshold);
    fprintf('\n');
    
    if (numel(S.Threshold)>1) && (numel(S.CellMask)==numel(S.Threshold))
        S.Threshold=median(S.Threshold(S.CellMask));
    end

    Str='    Final threshold applied:';
    fprintf('%-50s',Str);
    fprintf('th = %0.3f (CV = %0.2f%%)\n',S.Threshold,S.CVofData);
    clear('ScaleJ','DataCV');
    
    %% Photon Counting flow limits
    if S.Denoising>=1
        DataPhotonCounting=DataIntensity/S.Gain>=LowestPixelPhotonCount;
        MaximaMask=bsxfun(@times,MaximaMask,DataPhotonCounting);
        clear('DataPhotonCounting');
    end
    clear('DataIntensity');
    
    if ~isempty(S.PhotonCountingFlow)
        MaximaMask=bsxfun(@times,MaximaMask,S.PhotonCountingFlow>LowestPixelPhotonCount);
    end
    
    %% Applying cell mask.
    if S.ApplyGlobalCellmask
        fprintf('%-50s','    Applying cell mask:');
        MaximaMask=bsxfun(@times,MaximaMask,S.CellMask);
        fprintf('done\n');
    end

    %% Spark blob segmentation.
    if S.UseWatershedSegment
        fprintf('%-50s','    Sparks separation with watershed:');
        MaximaMask=CaEventSegmentWatershed(MaximaMask,S.xyt_dim)>0;
        fprintf('done\n');
    else
        if S.ApplyImdilation
            MaximaMask=imdilate(MaximaMask,ones(3,3,3));
            MaximaMask=imfill(MaximaMask,'holes');
            MaximaMask=imerode(MaximaMask,ones(3,3,3));
        end
    end

    fprintf('%-50s','    Sparks segmentation:');
        SparkLabel=bwconncomp(MaximaMask,ones(3,3,3));
    fprintf('done\n');
    
    % Applying spark volume limits here.
    fprintf('%-50s','    Applying spark volume limits:');
    area = cellfun(@numel, SparkLabel.PixelIdxList);
    DetectionLimit=S.DetectionLimit/S.xyt_dim(1)/S.xyt_dim(2)/S.xyt_dim(3);
    idxToKeep=(area>=DetectionLimit(1)) & (area<=DetectionLimit(2));
    SparkLabel.PixelIdxList=SparkLabel.PixelIdxList(idxToKeep);
    SparkLabel.NumObjects=sum(idxToKeep);
    fprintf('done\n');
    
    fprintf('%-50s','    Locating and labeling sparks:');
    S.SparkLabel=labelmatrix(SparkLabel);
    S.SparkArea=regionprops(S.SparkLabel,'Area','BoundingBox','Centroid');

    clear('SparkLabel','SparkArea','MaximaMask','DetectionLimit','area','idxToKeep');
    fprintf('%g potential sparks\n',numel(S.SparkArea));
    
    %% Information extraction of single sparks.
    fprintf('%-50s','    Single spark checking (info extraction):');
    [S.SparkPos,S.SparkProperty,S.SpatiotemporalFitting]=...
        SingleSparkv1_2_Parallel(S.Data,S.SparkArea,...
        S.SparkLabel,S.SparkLim,S.xyt_dim,S.CellMask);
    disp('done');

    %% Checking fitting results.
    fprintf('%-50s','    Check sparks with reasonable parameters:');
    [S.SparkPos,S.SparkProperty,S.SparkPosRejected,...
        S.SparkPropertyRejected]=CheckReasonableParameters(...
        S.SparkPos,S.SparkProperty,S.GauAmpLim,...
        S.TauLim,S.FWHMLimit);
    
    S.SpatiotemporalFittingRejected=S.SpatiotemporalFitting;
    S.SpatiotemporalFittingRejected(S.SparkPos(:,end),:)=[];
    S.SpatiotemporalFitting=S.SpatiotemporalFitting(S.SparkPos(:,end),:);
    
    fprintf('\n');

    clear('SparkMaskKept','CellMaskEdge','SparkMaskRejected');

    %% Spark site analysis.
    if size(S.SparkProperty,1)>=2
        Str='    Spark clustering:';
        fprintf('%-50sMin distance = %0.2fum',Str,S.MinSparkSiteDistance);
        [S.SparkSitePos,S.SparkSiteProperty,S.Spark_Site_Relation,~]=SparkSites(S.SparkPos,...
            S.SparkProperty,S.xyt_dim,'MinDistance',S.MinSparkSiteDistance);
        fprintf('\n');
    elseif size(S.SparkProperty,1)==1
        S.SparkSitePos=S.SparkPos;
        S.SparkSiteProperty=S.SparkProperty;
        S.Spark_Site_Relation.SparkSiteNo=1;
        S.Spark_Site_Relation.SparkID='1';
    elseif size(S.SparkProperty,1)<1
        S.SparkSitePos=S.SparkPos;
        S.SparkSiteProperty=S.SparkProperty;
        S.Spark_Site_Relation.SparkSiteNo=0;
        S.Spark_Site_Relation.SparkID='';
    end
    Str='    Spark site number:'; fprintf('%-50s%d final spark sites\n',Str,size(S.SparkSiteProperty,1));
    
    
    %% Spark/Site Frequency (in /squm/s).
    CellArea=sum(S.CellMask(:)*S.xyt_dim(1)*S.xyt_dim(2));
    TemporalDuration=RawSize(3)*S.xyt_dim(3)/1000;
    S.SparkFrequency=size(S.SparkPos,1)/TemporalDuration/CellArea;
    S.SparkSiteFrequency=size(S.SparkSitePos,1)/TemporalDuration/CellArea;
    clear('TemporalDuration','CellArea');
    
    %% Print elapsed time.
    Ctime=round(clock);
    Ctime=[sprintf('%02d',Ctime(4)) ':' sprintf('%02d',Ctime(5)) ':' sprintf('%02d',Ctime(6)) ...
        ', ' sprintf('%4d',Ctime(1)) '-' sprintf('%02d',Ctime(2)) '-' sprintf('%02d',Ctime(3))];
    fprintf('%-50s%s\n','    Analysis finished at:',Ctime);
    clear('Ctime')
    
    elaspsedTime=toc(totaltime)/60;
    if elaspsedTime>1
        Str='minutes';
    else
        Str='minute';
    end
    fprintf('%-50s%0.2f %s\n','    Total elapsed time:',elaspsedTime,Str);
    
    %% Clear unnecessary fields.
    S=rmfield(S,{'SparkArea'});
    
    %% Output
    if nargout==1
        varargout(1)={S};
    end
end
