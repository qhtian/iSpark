function varargout=SimCalAdpTH(varargin)
    % Simulate spark recordings and calculate suitable adaptive thresold.
    %
    % The simulation of spark recordings will be 500 frames with 120x512
    %   pixels. Pixel dimension is 0.28um x 0.28um x 10ms.
    %
    % Usage:
    %   FinalThNum=SimCalAdpTH('Key',value,...);
    %   The following lists all the keys and the default values:
    %     'Data',           [];  Raw recordings without Ca Spark to optimize speicific device.
    %     'Gain',           3 repeats of [0.2,0.4,0.7,round((1:0.5:20).^2)];
    %     'GaussSigma',     3.5;
    %     'Threshold',      (12:-0.5:1);
    %     'xyt_dim',        [0.28,0.28,10] in um,um,ms;
    %     'DetectionScale', [20,60] in ms;
    %     'DetectionLimit', [1,60000], Detection mass limit, e.g. 0.28um*0.28um*10ms=0.784.
    %     'RecordLength',   500, recording frames;
    %     'MaxFalsePositive',1500; Max final False Positive sparks.
    %
    %   All the results will be automatically saved as text file
    %   (SimCalAdpTH.txt), they can also be saved as variable.
    %
    % The field SimCalResults in FinalThNum has the following format:
    %   FinalThNum(1,:): input gain value;
    %   FinalThNum(2,:): rawCV;
    %   FinalThNum(3,:): detected gain value;
    %   FinalThNum(4,:): denoised CV;
    %   FinalThNum(5:end,:): threshold/potentialNum.
    %   FinalThNum(5:end,1:2:end): local threshold used;
    %   FinalThNum(5:end,2:2:end): potential number.
    %
    %   The FinalThNum can be analyzed by the AnalyzeAdaptiveThreshold.
    %
    % Example Usage:
    %   From real recordings:
    %     FinalThNum=SimCalAdpTH('Data',ImagestackFromRecording,'xyt_dim',[0.28,0.28,10]);
    %   Simulations of common PMT detecor with moderate settings:
    %     FinalThNum=SimCalAdpTH('Gain',150,'xyt_dim',[0.28,0.28,10]);
    %   Simulations of sCMOS camera, e.g. Hamamatsu Flash4.0:
    %     FinalThNum=SimCalAdpTH('Gain',2.5,'xyt_dim',[0.215,0.215,6.85]);
    %
    % See also AnalyzeAdaptiveThreshold.

    
    %% Input & parameters.
    % Gain=[0.2,0.4,0.7,round((1:0.5:14.5).^2)];
    Gain=200./[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2,3,5,8,11,15,20,30,50,100,200,300,500,1000];
    p=inputParser;
    % Experiment-specific simulation.
    p.addParameter('Gain',       Gain,@(x) sum((x<=1000) & (x>=0.2))>=1);  clear('Gain')
    p.addParameter('GaussSigma', 3.5,@(x) (x>0) && (x<=4095));
    p.addParameter('LaserIntensity', 40,@(x) (x>0) && (x<=4095));
    p.addParameter('xyt_dim',    [0.28,0.28,10],@(x)numel(x)==3 && min(x)>0);
    
    p.addParameter('RecordLength',500,@(x)numel(x)==1 && (x>0));
    
    % General settings.
    p.addParameter('DetectionScale', [20,10,60],@(x)numel(x)==3 && min(x)>0);
    p.addParameter('DetectionLimit',[6.27, 60000],@(x)numel(x)==2 && min(x)>0);%Detection mass limit.

    % Running options.
    p.addParameter('MaxFalsePositive',1500,@(x)numel(x)==1 && (x>0));
    p.addParameter('Threshold',  [25:-2:14,12:-0.5:1],@(x) (min(x)>0) && (max(x)<=10));

    % Cell Area, for sim without real cell
    p.addParameter('CellArea',   2472,@(x) numel(x)==1 && (x>5)); % um2.
    p.addParameter('BasalIntDen', 200,@(x) numel(x)==1 && (x>0));
    p.addParameter('Offset', 50,@(x) numel(x)==1 && (x>=0));
    
    % Switch.
    p.addParameter('SimWithRealCell',false,@(x) numel(x)==1 && islogical(x)); % use real cell or pure simulation.
    
    % Real recording without Spark as testing data.
    p.addParameter('Data', [],@(x)size(x,1)>80 && size(x,2)>80 && size(x,3)>300); % use with real Bgr recording.
    p.addParameter('CellMask', [],@(x)size(x,1)>80 && size(x,2)>80 && islogical(x)); % use with real Bgr recording.
    
    parse(p, varargin{:});
    p=p.Results;
    p.Threshold=    sort(p.Threshold,'descend');
    clear('varargin')
    
    %% Parameters for spark detection.
    Gain=p.Gain;
    Gain_num=numel(Gain);
    Gain=Gain((Gain>=0.2) & (Gain<=1000));
    if Gain_num>numel(Gain)
        disp('    Gain out of the range [0.2,1000] will be discarded.');
    end
    clear('Gain_num');
    
    GaussSigma=     p.GaussSigma;
    Threshold=    	p.Threshold;
    xyt_dim=        p.xyt_dim;
    DetectionScale= p.DetectionScale;
    DetectionLimits=p.DetectionLimit/xyt_dim(1)/xyt_dim(2)/xyt_dim(3);
    RecordLength=   p.RecordLength;
    MaxFalsePositive=p.MaxFalsePositive;
    Offset=         p.Offset;
    LaserIntensity=p.LaserIntensity;
    CellArea=       p.CellArea;
    BasalIntDen=    p.BasalIntDen;
    SimWithRealCell=p.SimWithRealCell;
    CellMask       =p.CellMask;
    %% Main cyles start here.
    if ~isempty(p.Data)
        Gain=NaN;
    end
    FinalThNum=zeros(numel(Threshold)+4,2*numel(Gain));
    for k=1:numel(Gain)
        currGainBlock=tic;
        % Progress indicator.
        disp(['    Gain: ' num2str(Gain(k)) ...
            '    No. ' num2str(k) ' of ' num2str(numel(Gain)) ...
            '    Progress: ' num2str(round(k/numel(Gain)*1000)/10) '%']);
        
        FinalThNum(1,(k*2-1):k*2)=Gain(k);
        
        % Spark recording simulaiton.
        if isempty(p.Data)
            if SimWithRealCell
                I=SparkRecordingSimWithRealBgr('Gain',Gain(k),'SparkNum',0,...
                    'Tdim',RecordLength,'SetupNoise',GaussSigma,'LaserIntensity',LaserIntensity,...
                    'SetupOffset',Offset);
                I=single(I);
                CellMask=MaskingIsodata(I);
            else
                x=sqrt(CellArea/4);
                y=x*4;
                
                x=round(x/xyt_dim(1));
                y=round(y/xyt_dim(2));
                I=zeros(x,y+x,RecordLength,'single');
                I(:,1:y,:)=I(:,1:y,:)+BasalIntDen;
                CellMask=false(x,y+x);
                CellMask(1:x,1:y)=true;
                
                parfor j=1:RecordLength
                    I(:,:,j)=Gain(k)*poissrnd(I(:,:,j)/Gain(k)); %#ok<PFBNS>
                end
                I=I+randn(size(I))*GaussSigma+Offset;
            end
        else
            I=single(p.Data);
        end
        
        % Cell Mask.
        % CellMask=CellMasking(I,1);
        % CellMask=imerode(CellMask,[0,1,0;1,1,1;0,1,0]);


        %% PURE denoise
        fprintf('%-50s','    PURE-LET Poissonian/Gaussian denoising:');
        DataCV=mad(diff(I,1,3),1,3)*1.4826;
        DataIntensity=mean(I,3);
        DataCV=DataCV./DataIntensity;
        DataCV(isnan(DataCV))=inf;
        DataCV=median(DataCV(CellMask));
        FinalThNum(2,(k*2-1):k*2)=DataCV;
        
        [I,DetectedGain,CameraOffset,~]=PUREDenoiseJava(I);
        I=I-CameraOffset;

     
        DataCV=mad(diff(I,1,3),1,3)*1.4826;
        DataIntensity=mean(I,3);
        DataCV=DataCV./DataIntensity;
        DataCV(isnan(DataCV))=inf;
        
        DataCV=median(DataCV(CellMask));
        FinalThNum(3,(k*2-1):k*2)=DetectedGain;

        FinalThNum(4,(k*2-1):k*2)=DataCV;
        
        fprintf('done\n');

        fprintf('\n\n    InputGain: %g    DetectedGain: %0.3f\tCV: %0.5f\n',Gain(k),DetectedGain,DataCV);
        clear('DetectedGain','FinalCV','ans','Str','s');
        
        %% Cycle Thresholds.
        for m=1:numel(Threshold)
            SmallBlock=tic;
            LocalThreshold=Threshold(m);
            fprintf('    Threshhold=%0.1f,',LocalThreshold);
            
            ScaleJ=round(DetectionScale/xyt_dim(3));
            % ScaleJ=[max(2,ScaleJ(1)),round(10/xyt_dim(3)),max(1,ScaleJ(2))];
            fprintf(' thresholding: ');
            [MaximaMask,~]=SparkThresholding(I,'ScaleN',ScaleJ,'Threshold',LocalThreshold);

            % Apply detection limits.
            fprintf('. Clean')
            MaximaMask=imclose(MaximaMask,strel('disk',1,0));
            MaximaMask=bsxfun(@and,MaximaMask,CellMask);
            
            % Segmentation.
            fprintf('. Segment')
            SparkLabel = bwconncomp(MaximaMask,ones(3,3,3));
            % Applying spark volume limits here.
            area = cellfun(@numel, SparkLabel.PixelIdxList);
            idxToKeep=(area>=DetectionLimits(1)) & (area<=DetectionLimits(2));
            PotentialNum=sum(idxToKeep);
            fprintf('. %g final sparks',PotentialNum);
            SmallBlock=toc(SmallBlock);
            disp(['. ' num2str(round(SmallBlock)) ' seconds.'])
            
            FinalThNum(4+m,k*2-1)=LocalThreshold;
            FinalThNum(4+m,k*2)=PotentialNum;
            clear('LocalThreshold','SparkLabel','SparkArea','MaximaMask','SparkArea','SparkPos',...
                'SparkProperty','ans');
            
            if PotentialNum>=MaxFalsePositive
                clear('Str','ans','I','FinalSparkNum')
                break
            end
            clear('PotentialNum')
        end

        %%
        clear('I_center','I_avg','I_std','LocalThreshold','PotentialNum','Str','ans','I')
        fprintf('  Current Gain Block: ');
        toc(currGainBlock)
        
        fprintf('\n\n\n')
    end
    
    %% Output
    p.SimCalResults=FinalThNum;
    if nargout==0
        assignin('base','SimCalAdpTH',p);
    end
    varargout(1)={p};
end


function CellMask=MaskingIsodata(I)
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













% function fitLinear(thresh,SpkNum)
%     SpkNum=log(SpkNum);
%     thresh_fit=zeros(size(SpkNum,2),1);
%     fprintf('   Exp-specific Threshold:  ')
%     for k=1:size(SpkNum,2)
%         currSpkNum=SpkNum(:,k);
%         currThresh=thresh(:,k);
%         bw=currSpkNum>0;
%         currSpkNum=currSpkNum(bw);
%         currThresh=currThresh(bw);
%         
%         p=robustfit(currThresh,currSpkNum);
%         thresh_fit(k)=-p(1)/p(2);
%         fprintf('\t%1.3f',thresh_fit(k));
%     end
%     fprintf('\n      Final Threshold:  %1.3f\n',median(thresh_fit));
% end


% function varargout=SimCalAdpTH(varargin)
%     % Simulate spark recordings and calculate suitable adaptive thresold.
%     %
%     % The simulation of spark recordings will be 500 frames with 120x512
%     %   pixels. Pixel dimension is 0.28um x 0.28um x 10ms.
%     %
%     % Usage:
%     %   FinalThNum=SimCalAdpTH('Key',value,...);
%     %   The following lists all the keys and the default values:
%     %     'LET_ID',         2;
%     %     'Gain',           3 repeats of [0.2,0.4,0.7,round((1:0.5:20).^2)];
%     %     'GaussDenoise',   'smoothn'; (option: 'smoothn' or 'Coif1');
%     %     'GaussSigma',     3.5;
%     %     'Threshold',      (5.6:-0.2:2.6)';
%     %     'xyt_dim',        [0.28,0.28,10] in um,um,ms;
%     %     'DetectionScale', [20,60] in ms;
%     %     'SparkLim',       [50 100]; (leading foot / tail length in ms);
%     %     'GauAmpLim',      10 in dFF0;
%     %     'TauLim',         500 in ms;
%     %     'FWHMLimit',      [0.2,10] in um, the Min/Max FWHM of spark.
%     %     'DetectionLimit', [6.27 60000], Detection mass limit.
%     %     'RecordLength',   500, recording frames;
%     %     'MaxFalsePositive',1500; Max final False Positive sparks.
%     %     'ExpSpecificThreshold', false; Simulate specific experiments. See
%     %                       later texts for details.
%     %
%     %   All the results will be automatically saved as text file
%     %   (SimCalAdpTH.txt), they can also be saved as variable.
%     %
%     %   Option: 'ExpSpecificThreshold'
%     %     When turn on this option, you can simulate and calculate experiment-speicifc
%     %     thereshold by yourself. For example:
%     %     SimCalAdpTH('Gain',YourDeviceGain,'GaussSigma',YourDeviceSigma,'ExpSpecificThreshold',true);
%     %     YourDeviceGain and YourDeviceSigma can be calculate by PURE_Denoise_Java function.
%     %     You can repeat time multiple times just by duplicating the input Gain.
%     %
%     % The field SimCalResults in FinalThNum has the following format:
%     %   FinalThNum(1,:): input gain value;
%     %   FinalThNum(2,:): rawCV;
%     %   FinalThNum(3,:): detected gain value;
%     %   FinalThNum(4,:): denoised CV;
%     %
%     %   FinalThNum(5:end,:): threshold/potentialNum.
%     %
%     %   FinalThNum(5:end,1:2:end): local threshold used;
%     %   FinalThNum(5:end,2:2:end): potential number.
%     %
%     %   The FinalThNum can be analyzed by the AnalyzeAdaptiveThreshold.
%     %
%     % See also AnalyzeAdaptiveThreshold.
%     %
%     % function Sample_Usage_Of_SimCalAdpTH
%     %     try
%     %         THGauss=SimCalAdpTH('GaussDenoise',0);    save('TH.mat','THGauss');
%     %         THnoGauss=SimCalAdpTH('GaussDenoise',-1);  save('TH.mat','THGauss','THnoGauss');
%     %     catch err
%     %         disp(err)
%     %         save('Stopped.mat');
%     %     end
%     %
%     %     eval('!shutdown -s -f -t 300') % Shutdown the PC in 5 minutes.
%     % end
%     %
%     % To shut down the computer:   eval('!shutdown -s -f -t 300');
%     % To cancel the shutting down: evalc('!shutdown -a');
%     
%     %% Input & parameters.
%     Gain=[0.2,0.4,0.7,round((1:0.5:14.5).^2)];
%     % Gain=cat(2,Gain,Gain,Gain)';
%     
%     p=inputParser;
%     
%     % Experiment-specific simulation.
%     p.addParameter('Gain',       Gain,@(x) sum((x<=1000) & (x>=0.2))>=1);  clear('Gain')
%     p.addParameter('GaussSigma', 3.5,@(x) (x>0) && (x<=4095));
%     p.addParameter('LaserIntensity', 40,@(x) (x>0) && (x<=4095));
%     p.addParameter('Threshold',  (8:-0.2:0.2)',@(x) (min(x)>0) && (max(x)<=10));
%     p.addParameter('xyt_dim',    [0.28,0.28,10],@(x)numel(x)==3 && min(x)>0);
%     p.addParameter('RecordLength',500,@(x)numel(x)==1 && (x>0));
%     p.addParameter('ExpSpecificThreshold',false,@(x)numel(x)==1 && islogical(x));
%     
%     % General settings.
%     p.addParameter('DetectionScale', [20,60],@(x)numel(x)==2 && min(x)>0);
%     p.addParameter('SparkLim',   [50,100],@(x)numel(x)==2 && min(x)>0); %[leading foot/tail length]
%     p.addParameter('GauAmpLim',  10,@(x)isscalar(x) && x>0);
%     p.addParameter('TauLim',     500,@(x)isscalar(x) && x>0);
%     p.addParameter('FWHMLimit',  [0.2,10],@(x)numel(x)==2 && min(x)>0);% Min/Max FWHM of spark.
%     p.addParameter('DetectionLimit',[1.00, 60000],@(x)numel(x)==2 && min(x)>0);%Detection mass limit.
%     
%     % Running options.
%     p.addParameter('MaxFalsePositive',1500,@(x)numel(x)==1 && (x>0));
%     p.addParameter('PureEmptyImageStack',false,@(x)numel(x)==1 && islogical(x));
%     
%     parse(p, varargin{:});
%     p=p.Results;
%     p.Threshold=    sort(p.Threshold,'descend');
%     clear('varargin')
%     
%     %% Parameters for spark detection.
%     Gain=p.Gain;
%     Gain_num=numel(Gain);
%     Gain=Gain((Gain>=0.2) & (Gain<=1000));
%     if Gain_num>numel(Gain);disp('    Gain out of the range [0.2,1000] will be discarded.');end
%     clear('Gain_num');
%     
%     GaussSigma=     p.GaussSigma;
%     Threshold=    	p.Threshold;
%     xyt_dim=        p.xyt_dim;
%     DetectionScale= p.DetectionScale;
%     % SparkLim=       p.SparkLim;
%     % GauAmpLim=      p.GauAmpLim;
%     % TauLim=         p.TauLim;
%     % FWHMLimit=      p.FWHMLimit;
%     DetectionLimits=p.DetectionLimit;
%     DetectionLimits=DetectionLimits/xyt_dim(1)/xyt_dim(2)/xyt_dim(3);
%     if DetectionLimits(1)<2; DetectionLimits(1)=2; end                % At least remove single pixels.
%     RecordLength=   p.RecordLength;
%     MaxFalsePositive=p.MaxFalsePositive;
%     PureEmptyImageStack=p.PureEmptyImageStack;
%     ExpSpecificThreshold=p.ExpSpecificThreshold;
%     LaserIntensity=p.LaserIntensity;
% 
%     %% Main cyles start here.
%     FinalThNum=zeros(numel(Threshold)+4,2*numel(Gain));
%     % logfilename='SimCalAdpTH.txt';
%     for k=1:numel(Gain)
%         currGainBlock=tic;
%         % Progress indicator.
%         disp(['    Gain: ' num2str(Gain(k)) ...
%             '    No. ' num2str(k) ' of ' num2str(numel(Gain)) ...
%             '    Progress: ' num2str(round(k/numel(Gain)*1000)/10) '%']);
%         
%         FinalThNum(1,(k*2-1):k*2)=Gain(k);
%         
%         if PureEmptyImageStack
%             I=poissrnd(200+zeros(120,512,RecordLength))+randn(120,512,RecordLength)*GaussSigma;%randn(120,512,RecordLength)*GaussSigma;
%             CellMask=true(120,512);
%         else
%             % Spark recording simulaiton.
%             I=SparkRecordingSimWithRealBgr('Gain',Gain(k),'SparkNum',0,...
%                 'Tdim',RecordLength,'SetupNoise',GaussSigma,'LaserIntensity',LaserIntensity);
%             I=double(I);
%             
%             % Cell Mask.
%             CellMask=CellMasking(I,1);
%             CellMask=imerode(CellMask,[0,1,0;1,1,1;0,1,0]);
%         end
% 
%         
%         
%         
%         %% PURE denoise
%         fprintf('%-50s','    PURE-LET Poissonian/Gaussian denoising:');
%         DataCV=mad(diff(I,1,3),1,3)*1.4826;
%         DataIntensity=mean(I,3);
%         DataCV=DataCV./DataIntensity;
%         DataCV(isnan(DataCV))=inf;
%         DataCV=median(DataCV(CellMask));
%         FinalThNum(2,(k*2-1):k*2)=DataCV;
%         
%         [I,DetectedGain,CameraOffset,~]=PUREDenoiseJava(I);
%         if ~PureEmptyImageStack
%             I=I-CameraOffset;
%         end
%      
%         DataCV=mad(diff(I,1,3),1,3)*1.4826;
%         DataIntensity=mean(I,3);
%         DataCV=DataCV./DataIntensity;
%         DataCV(isnan(DataCV))=inf;
%         
%         DataCV=median(DataCV(CellMask));
%         FinalThNum(3,(k*2-1):k*2)=DetectedGain;
% 
%         FinalThNum(4,(k*2-1):k*2)=DataCV;
%         
%         fprintf('done\n');
%         %% Write Gain and CV to file.
%         % filelog = fopen(logfilename,'a+');
%         % fprintf(filelog,'InputGain: %g\tDetectedGain: %0.3f\tCV: %0.5f',Gain(k),DetectedGain,DataCV);
%         % fclose(filelog);
%         fprintf('\n\n    InputGain: %g    DetectedGain: %0.3f\tCV: %0.5f\n',Gain(k),DetectedGain,DataCV);
%         clear('DetectedGain','FinalCV','ans','Str','s');
%         
%         %% Cycle Thresholds.
%         for m=1:numel(Threshold)
%             SmallBlock=tic;
%             LocalThreshold=Threshold(m);
%             fprintf('    Threshhold=%0.1f,',LocalThreshold);
%             
%             ScaleJ=round(DetectionScale/xyt_dim(3));
%             ScaleJ=[max(2,ScaleJ(1)),round(10/xyt_dim(3)),max(1,ScaleJ(2))];
%             fprintf(' thresholding: ');
%             [MaximaMask,~]=SparkThresholding(I,'ScaleN',ScaleJ,'Threshold',LocalThreshold);
% 
%             % Apply detection limits.
%             fprintf('. Clean')
%             MaximaMask=imclose(MaximaMask,strel('disk',1,0));
%             MaximaMask=bsxfun(@and,MaximaMask,CellMask);
%             
%             % Segmentation.
%             fprintf('. Segment')
%             % MaximaMask=MaximaMask.*repmat(CellMask,[1 1 size(I,3)]);
%             SparkLabel = bwconncomp(MaximaMask,ones(3,3,3));
%             % Applying spark volume limits here.
%             area = cellfun(@numel, SparkLabel.PixelIdxList);
%             idxToKeep=(area>=DetectionLimits(1)) & (area<=DetectionLimits(2));
%             % SparkLabel.PixelIdxList=SparkLabel.PixelIdxList(idxToKeep);
%             % SparkLabel.NumObjects=sum(idxToKeep); clear('idxToKeep')
%             
%             % SparkLabel = labelmatrix(SparkLabel);
%             % SparkArea = regionprops(SparkLabel,'Area','BoundingBox','Centroid');
% 
%             %Spark Information extraction of single sparks.
%             % fprintf('. Extract')
%             % [SparkPos,SparkProperty]=SingleSparkv1_2_Parallel(I,SparkArea,SparkLabel,SparkLim,xyt_dim,CellMask);
%             
%             % Checking fitting results.
%             % fprintf('. Check ')
%             % [SparkPos,~,~,~,]=CheckReasonableParameters(SparkPos,SparkProperty,GauAmpLim,TauLim,FWHMLimit);
%             
%             PotentialNum=sum(idxToKeep);
%             fprintf('. %g final sparks',PotentialNum);
%             SmallBlock=toc(SmallBlock);
%             disp(['. ' num2str(round(SmallBlock)) ' seconds.'])
% 
%             % PotentialNum=size(SparkPos,1);
%             
%             FinalThNum(4+m,k*2-1)=LocalThreshold;
%             FinalThNum(4+m,k*2)=PotentialNum;
%             % filelog = fopen(logfilename,'a+');
%             % fprintf(filelog,'\tThreshold: %0.1f\tFinalSparkNum: %g',LocalThreshold,PotentialNum);
%             % fclose(filelog);
%             clear('LocalThreshold','SparkLabel','SparkArea','MaximaMask','SparkArea','SparkPos',...
%                 'SparkProperty','ans');
%             
%             if PotentialNum>=MaxFalsePositive
%                 clear('Str','ans','I','FinalSparkNum')
%                 break
%             end
%             clear('PotentialNum')
%         end
%         % Cycle threshold finished here.
%         %%
%         clear('I_center','I_avg','I_std','LocalThreshold','PotentialNum','Str','ans','I')
%         % filelog = fopen(logfilename,'a+');
%         % fprintf(filelog,'\n');
%         % fclose(filelog);
%         fprintf('  Current Gain Block: ');
%         toc(currGainBlock)
%         
%         fprintf('\n\n\n')
%     end
%     
%     %% Output
%     if nargout==1
%         p.SimCalResults=FinalThNum;
%         varargout(1)={p};
%     end
%     
%     
%     %% ExpSpecificThreshold.
%     if ExpSpecificThreshold
%         fitLinear(FinalThNum(5:end,1:2:end),FinalThNum(5:end,2:2:end));
%     end
% end
% 
% 
% function fitLinear(thresh,SpkNum)
%     SpkNum=log(SpkNum);
%     thresh_fit=zeros(size(SpkNum,2),1);
%     fprintf('   Exp-specific Threshold:  ')
%     for k=1:size(SpkNum,2)
%         currSpkNum=SpkNum(:,k);
%         currThresh=thresh(:,k);
%         bw=currSpkNum>0;
%         currSpkNum=currSpkNum(bw);
%         currThresh=currThresh(bw);
%         
%         p=robustfit(currThresh,currSpkNum);
%         thresh_fit(k)=-p(1)/p(2);
%         fprintf('\t%1.3f',thresh_fit(k));
%     end
%     fprintf('\n      Final Threshold:  %1.3f\n',median(thresh_fit));
% end