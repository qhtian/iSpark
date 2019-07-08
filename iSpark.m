function iSpark(xyt_dim)
    % iSpark: A graphic interface for the analysis of cardiac calcium release
    % events, the calcium sparks.
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
    % Gaussian noise Handling: discrete consine transform-based GCV minimization or
    % zeroing fisrt level of wavelet (Coif1) decomposition. For minizaing GCV score,
    % See paper from Garcia D, Robust smoothing of gridded data in one and higher
    % dimensions with missing values. Computational Statistics & Data Analysis, 2010.
    %
    %
    % If the intermediate result is saved, it has the following fields:
    %  1.  DetectionLimit:     the same parameter from the input.
    %  2.  GauAmpLim:   the same parameter from the input.
    %  3.  Threshold:   the same parameter from the input.
    %  4.  SparkLim:    the same parameter from the input.
    %  5.  xyt_dim:     the same parameter from the input.
    %  6.  Creator:     the program used to generate the results.
    %  7.  Data:        the denoised data of input image stack.
    %  8.  Gain:        the original device gain value.
    %  9.  CellMask:    detected cell mask.
    %  10. CellMaskEdge: detected cell edge from CellMask.
    %  11. CameraOffset: dark background of the detector.
    %  12. GaussNoise:   sigma of dark noise.
    %  13. SparkLabel:   all the sparks detected and labeled in specific number.
    %  14. SparkPos: the coordinates of the final sparks, in the format of
    %      [x1,x2,y1,y2,t1,t2,xc,yc,t_onset,ID_In_SparkLabel].
    %  15. SparkProperty: the calculated spark properties, in the format of
    %      [FWHM,FWHM_R2,Amp,Bgr,Tau,Sigma,Amp_R2,Old_ID];
    %  16. SparkSitePos: same information as SparkPos, but for spark sites.
    %  17. SparkSiteProperty: same information as SparkProperty, but for spark sites.
    %  18. Spark_Site_Relation: the mapping of single sparks to spark sites.
    %  19. MaximaMask_1_Raw: Detected spark mask in first round.
    %  20. SparkPosRejected: the coordinates of the rejected sparks, in the format
    %      of [x1,x2,y1,y2,t1,t2,xc,yc,t_onset,ID_In_SparkLabel].
    %  21. SparkPropertyRejected: the rejected spark properties, in the format of
    %      [FWHM,FWHM_R2,Amp,Bgr,Tau,Sigma,Amp_R2,Old_ID,AmpCheck,FWHMCheck,MassCheck];
    %  22. MaximaMask_2_CheckingKept: finally detected spark masks.
    %  23. MaximaMask_2_CheckingRejected: detected but finally rejected spark masks.
    %  24. SparkFrequency: spark frequency and spark site frequency (/squm/s).
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
    %
    %
    % Change log v1.2:
    %   1. Live view switch.
    %   2. Hard-wired numbers are changed to relative numbers based on xyt_dim.
    %   3. Adjusted live view.
    %   4. Peak frame fitting is now with Intensity data, not F/F0;
    %   5. All the outputs are now in um, squm or ms.
    % Change log v1.32:
    %   1. Reconstruction is now based on a all-in-one 3D spark formula.
    %   2. A spark movie simulation program is now available.
    % Change log v1.33:
    %   1. Noise estimation is now performed with the model of gain.
    % Change log v1.34:
    %   1. Preliminary adaptive threholding to remove poissoninan noise.
    %   2. Spark simulation removed.
    % Change log v1.35:
    %   1. Refined adaptive threholding to remove poissoninan noise.
    %   2. Refined spark simulation with GauConvExp.
    %   3. Rewritten codes for the extraction of single spark inforation.
    % Change log v1.36:
    %   1. New function module: result checking.
    % Change log v1.37:
    %   1. Validated 'LocalThreshold' values for all kinds of gain.
    % Change log v1.38:
    %   1. Refined single spark info extraction.
    %   2. Spark mapping function.
    % Change log v1.39:
    %   1. Adaptive local thresholding.
    % Change log v1.40:
    %   1. Refined cell mask detection.
    % Change log v1.41:
    %   1. Refined Gain determination and Gain-Threshold relationship.
    %   2. Refined spark mass calculation with all-in-one formula.
    % Change log v1.42:
    %   1. Last hard-wired numbers removed.
    %   2. Local single pixel number is now changed to relative value.
    % Change log v1.43 BranchA:
    %   1. Adaptive Gaussian denoising, which make it suitable for both
    %      HyD data, PMT data and camera data.
    %   2. Threshold is based on the coefficient of variation after data
    %      denoising.
    % Change log v1.45:
    %   1. Threshold is determined for each single pixel.
    %   2. Completely adaptive denosing for both poissnian noise and
    %      Gaussain noise (basal Ca vibrations).
    % Change log v1.46:
    %   1. Advanced mode in the SparkFolderAnalysis.
    %   2. New user interface: SparkAnalysis.
    %
    % Change log v2.20:
    %   1. PUREshrink denoising.
    %
    % Change log v3.00:
    %   1. PURE-LET denoising.
    %
    % Change log v4.00:
    %   1. Spark thresholding using A Trous wavelet tranform-like upstroke detection.
    %   2. Improved processing speed.
    % Change log v4.10:
    %   1. Sensitivity choices.
    % Change log v4.20:
    %   1. Parallel computation supports.
    % Change log v4.30:
    %   1. Complete parallel computation, with reduced memory requirement.
    % Change log v4.40:
    %   1. Organized codes and anotations.
    %   2. Optimized isodata algorithm.
    %   3. Debugging for Zeiss LSM 5 LIVE data.
    % Change log v4.50:
    %   1. Adaptive Gaussian noise removal with smoothn.
    % Change log v4.60:
    %   1. Supports with bio-formats.
    % Change log v4.70:
    %   1. Manual mask checking/drawing complementary to automatic detection.
    %   2. Manual detector offset setting.
    %   3. Combined output files for spark information.
    %   4. Reorganized batch processing with the support of subfolders.
    %   5. Optimized calculations of offset and noise variance.
    % Change log v4.80:
    %   1. Phase contrast shift correction for point-laser scanning microscopy.
    % Change log v4.90:
    %   1. Change the code of PURE-LET denoising.
    %   2. Change the code for linear fitting of the device gain.
    % Change log v4.97:
    %   1. Change of smoothn to minimize |(Ioriginal-Idenoised)^2-MSE|, which is more reliable.
    %   2. Change the code of DeviceGain to get more realistic and faster results.
    %   3. Cover the data from sCMOS camera.
    %   4. Cover the data of neuron cells.
    %   5. Cover the data of nanosparks.
    % Change log v4.98:
    %   Quite a lot small changes.
    % Change log v4.99:
    %   1. Output changes.
    % Change log v5.00:
    %   1. Change the denoise to PURE-LET Java Plugin for ImageJ.
    % Change log v5.01:
    %   1. Support multiple cells.
    % Change log v5.02:
    %   1. Watershed to split connected multiple events.
    %   2. New image denoise method CANDLEdenoise added.
    %
    % Developed by Qinghai Tian <Qinghai.Tian@uniklinikum-saarland.de>.
    %      Institute for Molecular Cell Biology and Research Center for Molecular
    %      Imaging and Screening, Medical Faculty, Building 61, Saarland University,
    %      Homburg/Saar, Germany 66421.
    %
    % For research use only. For internal use only.
    
    
    
    %% Main Interface.
    [S.figSize,H,W,S.Screen]=FigSize;
    S.fh = figure('units','normalized','position',S.figSize,'menubar','none',...
        'name','iSpark - Automatic and intelligent calcium spark analysis',...
        'numbertitle','off','resize','on','Color',[0.85 0.85 0.85],...
        'HandleVisibility','off');
    
    %% Caption - Title.
    S.tx = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[60/W 638/H 400/W 30/H],'string','iSpark','fontweight','bold',...
        'backgroundcolor',get(S.fh,'color'),'fontsize',18,'FontUnits','Normalized');
    
    %% About button.
    S.pbAbout = uicontrol('Parent',S.fh,'style','push','unit','normalized',...
        'position',[460/W 620/H 30/W 17/H],'FontUnits','Normalized',...
        'string','Help','callback',{@About,S},...
        'tooltip','See features and operation precedures.');
    
    %% Caption - Basic Information.
    S.txFolder = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[30/W 593/H 200/W 25/H],'fontsize',13,'FontUnits','Normalized',...
        'HorizontalAlignment','left','fontweight','bold',...
        'string','Basic Information','backgroundcolor',get(S.fh,'color'));
    
    %% Folders to analyze.
    S.txFiles = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 570/H 300/W 25/H],'string','Files/Folders to Analyze:',...
        'backgroundcolor',get(S.fh,'color'),'fontsize',10,'FontUnits','Normalized',...
        'HorizontalAlignment','left');
    S.ls = uicontrol('Parent',S.fh,'style','list','unit','normalized',...
        'position',[50/W 420/H 310/W 156/H],'fontsize',10,'FontUnits','Normalized',...
        'string',{},'BackgroundColor',[1,1,1],'tooltip',sprintf(['Items that will be analyzed.\n \n' ...
        'Items can be added or deleted via\nthe buttons on the right side.']));
    
    S.pbAddFile = uicontrol('Parent',S.fh,'style','push','unit','normalized',...
        'position',[380/W 552/H 50/W 25/H],'FontUnits','Normalized',...
        'string','Add','callback',{@AddTiffFile,S},'fontWeight','bold',...
        'tooltip',sprintf(['Click to add items to the list of the ' ...
        'left listbox for analysis:\n \n1. *.tif or *.tiff image stacks. Or,' ...
        '\n2. folders containing *.tif or *.tiff images stacks.']));
    S.pbDelFolder = uicontrol('Parent',S.fh,'style','push','unit','normalized',...
        'position',[440/W 552/H 50/W 25/H],'FontUnits','Normalized',...
        'string','Delete','callback',{@DeleteItem,S},'fontWeight','bold',...
        'tooltip','Click to delete the item selected in the left listbox.');
    
    %% Save button.
    S.txSave = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 380/H 50/W 25/H],...
        'string','Output:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',10,'FontUnits','Normalized','HorizontalAlignment','left');
    S.Save = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[50/W 370/H 310/W 18/H],'backgroundcolor',[1,1,1],'enable','off',...
        'fontsize',9,'FontUnits','Normalized','string',...
        'Click to change result saving folder...',...
        'tooltip','Set the output file for the final results.','HorizontalAlignment','left');
    set(S.Save,'ButtonDownFcn',{@CreateSaveFile,S})
    
    
    %% xyt_dim.
    if nargin==0
        xyt_dim=[0.28,0.28,10];
    end
    S.txSpatialDim = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[380/W 515/H 250/W 25/H],'fontsize',10,...
        'string','Tempospatial Info:','backgroundcolor',get(S.fh,'color'),...
        'FontUnits','Normalized','HorizontalAlignment','left');
    S.xyt_dim1 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[380/W 507/H 32/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',num2str(xyt_dim(1)),'backgroundcolor',[1,1,1],...
        'tooltip','Set the spatial dimension (in um) for each pixel.');
    S.xyt_dim2 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[418/W 507/H 32/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',num2str(xyt_dim(2)),'backgroundcolor',[1,1,1],...
        'tooltip','Set the spatial dimension (in um) for each pixel.');
    S.xyt_dim3 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[456/W 507/H 32/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',num2str(xyt_dim(3)),'backgroundcolor',[1,1,1],...
        'tooltip','Set the temporal dimesion (in ms) for each frame.');

    %% Cell Masking algorithms.
    S.txCellMasking = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[380/W 468/H 150/W 25/H],'fontsize',10,...
        'string','Cell Masking:','backgroundcolor',get(S.fh,'color'),...
        'FontUnits','Normalized','HorizontalAlignment','left');
    S.CellMasking = uicontrol('Parent',S.fh,'style','pop','unit','normalized',...
        'position',[380/W 460/H 110/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',{'Isodata Algorithm';'Local Adaptive Thresholding';'Mask >= 3 * std ';...
        'Skip'},'backgroundcolor',[1,1,1],'value',1,...
        'tooltip',sprintf(['Cell Masking methods.\n' ...
        '--------------------------------------\n' ...
        '\nIsodata: Good for cardiomyocttes.',...
        '\nLocal Adapt.: Good for neuron cells.',...
        '\nMask >= 3*std: signal out of backgournd noise.',...
        '\nSkip: Do not try to mask cells.']));
    

    % Manual ROI drawing
    S.chROIDrawing = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [380/W 425/H 120/W 20/H],'string','Manual Mask Check','value',1,'FontUnits','Normalized','backgroundcolor',...
        get(S.fh,'color'),'tooltip',sprintf(['Check this on to manually check or draw cell area,\n' ...
        'e.g. in the case of permeabilized cell recordings.']));
    
    S.chWatershed = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [380/W 405/H 120/W 20/H],'string','Use Watershed','value',1,'FontUnits','Normalized','backgroundcolor',...
        get(S.fh,'color'),'tooltip','Use watershed algorithm to segment connected events.');
    
    S.chROIApplying = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [380/W 385/H 120/W 20/H],'string','Use Mask in Detect','FontUnits',...
        'Normalized','backgroundcolor',get(S.fh,'color'),'tooltip',...
        sprintf(['Check this on to apply detected cell mask to the\n' ...
        'detection area to remove false positives in background area.\n',...'
        'False positives could be significant in PMT/HyD data.']),'value',1);
    
    S.chImdilationApplying = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [380/W 365/H 120/W 20/H],'string','Use imdilate in Detect','FontUnits',...
        'Normalized','backgroundcolor',get(S.fh,'color'),'tooltip',...
        sprintf('Check this on to apply imdilate/imerode on detected spark masks'),'value',0);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%       Advanced mode.      %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S.txAdvancedMode = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [30/W 326/H 470/W 25/H],'string','Parameters for Advanced Mode','backgroundcolor',get(S.fh,'color'),...
        'fontsize',13,'FontUnits','Normalized','HorizontalAlignment','left','fontweight','bold');
    
    %% PURE-LET denoising ID.
    S.txImageDenoise = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 295/H 150/W 25/H],...
        'string','Image Denoise:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.ImageDenoise = uicontrol('Parent',S.fh,'style','pop','unit','normalized',...
        'position',[167/W 306/H 100/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',{'PURE-LET Denoise';'CANDLE Denoise';'Skip'},'backgroundcolor',...
        get(S.fh,'color'),'value',1,...
        'tooltip',sprintf(['PURE-LET denoising to remove Poissonian and Gaussian noise.\n' ...
        '--------------------------------------------------------\n' ...
        '\nPerform: perform the denoising process.',...
        '\nSkip: skip this process if the data is already denoised.']));
    
    %% DetectionScale.
    S.txDetectionScale = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [50/W 266/H 150/W 25/H],'string','Detect Scale (ms):','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.DetectionScaleMin = uicontrol('Parent',S.fh,'style','edit','unit','normalized','position',...
        [167/W 277/H 45/W 17/H],'string','20','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized',...
        'tooltip','Set roughtly to the temporal resolution or larger number.');
    S.DetectionScaleMax = uicontrol('Parent',S.fh,'style','edit','unit','normalized','position',...
        [222/W 277/H 45/W 17/H],'string','60','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized',...
        'tooltip','Set roughtly to two times of typical spark upstroke duration.');
    
    %% Theshold button.
    S.txLocalThreshold = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [295/W 295/H 150/W 25/H],'string','Spark Threshold:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.LocalThreshold = uicontrol('Parent',S.fh,'style','edit','unit','normalized','position',...
        [400/W 304/H 87/W 17/H],'string','','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized',...
        'tooltip',sprintf(['Set the threshold for local spark masking.\n \n' ...
        'Any specific threshold can be defined (recommended in 2~6).']));
    
    %% Spark site analysis parameter.
    S.txSparkSiteMinDistance = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [295/W 266/H 150/W 25/H],'string','Min. Site Distance:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.SparkSiteMinDistance = uicontrol('Parent',S.fh,'style','edit','unit','normalized','position',...
        [400/W 277/H 87/W 17/H],'string','1.0','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','tooltip',...
        sprintf(['Set the minimum distance (in um) that two sparks are\n'...
        'considered to originate from different RyRs.']));

    
    
    %% Dashed line to seperate different advanced parameters.
    annotation(S.fh,'line',[50/W 488/W],[265/H 265/H],'LineStyle',':');
    
    
    %% Amplitude limit button.
    S.txAmpLim = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 225/H 150/W 25/H],'fontsize',9,'FontUnits','Normalized',...
        'string','Max Amp Limit (dF/F0):','backgroundcolor',get(S.fh,'color'),...
        'HorizontalAlignment','left');
    S.AmpLim = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[170/W 233/H 106/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','10.0','backgroundcolor',get(S.fh,'color'),...
        'tooltip',sprintf('Set the maximum amplitude (in dF/F0)\nthat individual sparks can have.'));
    
    %% Decay tau limit.
    S.txDecayTau = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[295/W 225/H 150/W 25/H],...
        'string',sprintf('Max. Decay Limit:'),'backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.DecayTau = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[400/W 233/H 87/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','1500','backgroundcolor',get(S.fh,'color'),'tooltip',...
        sprintf('Set the maximum decay constant (the\ntau, in ms) that individual sparks can have.'));
    
    
    %% Offset
    S.txOffset = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[295/W 200/H 150/W 25/H],...
        'string',sprintf('Detector Offset:'),'backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.Offset = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[400/W 208/H 87/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','','backgroundcolor',get(S.fh,'color'),'tooltip',...
        sprintf('Set the detector / camera offset here.'));
    
    
    %% RecordPhase
    S.txRecordPhase = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[295/W 175/H 150/W 25/H],...
        'string',sprintf('Pahse Length:'),'backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.RecordPhase = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[400/W 183/H 87/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','','backgroundcolor',get(S.fh,'color'),'tooltip',...
        sprintf('Set the detector / camera offset here.'));

    %% Detection Limits.
    S.txDetectionLim = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [50/W 199/H 200/W 25/H],'string','Spark Volume Limits:','backgroundcolor',...
        get(S.fh,'color'),'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.DetectionLim1 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'fontsize',9,'FontUnits','Normalized','position',...
        [170/W 207/H 51/W 17/H],'string','6.27','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the minimum detection volume (um x um x ms) for single spark.');
    S.DetectionLim2 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'fontsize',9,'FontUnits','Normalized','position',...
        [225/W 207/H 51/W 17/H],'string','60000','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the maximum detection volume (um x um x ms) for single spark.');

    %% FWHM Limits.
    S.txFWHMLim = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 174/H 200/W 25/H],...
        'string','FWHM Limits (um):','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.FWHMLim1 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[170/W 182/H 51/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','0.2','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the minimum Full Width at Half Maximum (in um) for single spark.');
    S.FWHMLim2 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[225/W 182/H 51/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','10.0','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the maximum Full Width at Half Maximum (in um) for single spark.');
    
    %% Leading / Tail limits.
    S.txTail = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 147/H 200/W 25/H],...
        'string','Leading / Tail (ms):','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.Tail1 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[170/W 155/H 51/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','50','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the baseline length (in ms) before the spark.');
    S.Tail2 = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[225/W 155/H 51/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','150','backgroundcolor',get(S.fh,'color'),'tooltip',...
        'Set the tail length (in ms) after the spark.');
    
    %% Phase Contrast
    S.txPhase = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 120/H 150/W 25/H],...
        'string','Phase Contrast:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.Phase = uicontrol('Parent',S.fh,'style','pop','unit','normalized',...
        'position',[170/W 131/H 106/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',{'No correction';'Along X Axis';'Along Y Axis'},'backgroundcolor',...
        get(S.fh,'color'),'value',1,'tooltip',...
        sprintf(['Only valid for point-laser scanning microscopy.\n' ...
        'This option is to correct interlaced shift during bidirectional recordings.\n' ...
        ' \nIf not necessary, please always try to avoid it.\n' ...
        '------------------------------------------------------------------\n' ...
        '\nNone: Do not correct.',...
        '\nX Axis: Correct along the X axis.',...
        '\nY Axis: Correct along the Y axis.']));
    
    %% Correction of Image translation.
    S.txImRegister = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 92/H 150/W 25/H],...
        'string','Image translation:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.ImRegister = uicontrol('Parent',S.fh,'style','pop','unit','normalized',...
        'position',[170/W 103/H 106/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string',{'No correction';'Correct it!'},'backgroundcolor',...
        get(S.fh,'color'),'value',1,'tooltip',...
        sprintf(['Correct image difting during long time recording.\n' ...
        'This option is to correct image shift during a time lapse recording.\n' ...
        '------------------------------------------------------------------\n' ...
        'If not necessary, please always try to avoid it, since:\n' ...
        '1. The calculation takes very long time;\n2. Better to avoid it physically.']));
    
    %% sCMOS dynamic offset correction.
    S.txsCMOSOffset = uicontrol('Parent',S.fh,'style','text','unit','normalized',...
        'position',[50/W 65/H 150/W 25/H],...
        'string','sCMOS Offset Correct:','backgroundcolor',get(S.fh,'color'),...
        'fontsize',9,'FontUnits','Normalized','HorizontalAlignment','left');
    S.sCMOSdynamicOffsetCorrection = uicontrol('Parent',S.fh,'style','edit','unit','normalized',...
        'position',[170/W 73/H 106/W 17/H],'fontsize',9,'FontUnits','Normalized',...
        'string','NaN, 0, 5, 95, 100','backgroundcolor',get(S.fh,'color'),...
        'tooltip',sprintf(['Set sCMOS dynamic offset corrections.\n',...
        'The right format should be: [dim, x1, x2, x3, x4].\n',...
        'The first element decides which dim should be used. NaN means skip. Only dim = 1 / 2 is possible.\n',...
        'x1:x2 dedines the first reference area, and x3:x4 for the second.\n',...
        'Note that x1, x2 x3 and x4 are in percentage.']));
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%       Debug and testing modes.      %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S.txRunningOption = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [30/W 35/H 150/W 25/H],'string','Running Options','backgroundcolor',get(S.fh,'color'),...
        'fontsize',13,'FontUnits','Normalized','HorizontalAlignment','left','fontweight','bold');
    %% Save all results button.
    S.chSaveAll = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [50/W 15/H 70/W 20/H],'string','Save ALL','fontsize',9,'FontUnits','Normalized',...
        'backgroundcolor',get(S.fh,'color'),...
        'tooltip',sprintf(['Check this on to save all the intermediate processing results.'...
        '\nSpark info is always saved and not affected by this option.' ...
        '\n \nAll the results will be saved to the same folder as the input file,' ...
        '\nand will not be affected by the path for spark information.']),'value',0);
    %% Shutdown the system after analysis.
    S.chShutdown = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [130/W 15/H 80/W 20/H],'string','Power OFF','fontsize',9,'FontUnits','Normalized',...
        'backgroundcolor',get(S.fh,'color'),...
        'tooltip',sprintf(['Shut down the system after all the analysis are finished.' ...
        '\n \nFor Windows system, it works fine without the super-user privilege.\n \n'...
        'For Mac OS X and Linux system, Matlab should be executed by the super-user\n' ...
        'who is supposed to have the privilege to shut down the system.\n \n'...
        'For Mac OS X or Linux system, in a terminal use ''sudo'' to start matlab:\n'...
        'sudo / (your-matlab-folder) / bin / matlab']));
    %% Cut the original image stack for faster processing.
    S.chIcrop = uicontrol('Parent',S.fh,'style','checkbox','unit','normalized','position',...
        [220/W 15/H 80/W 20/H],'string','Corp Image','fontsize',9,'FontUnits','Normalized',...
        'backgroundcolor',get(S.fh,'color'),'Value',0,...
        'tooltip',sprintf('Corp the image for fast processing if a mask is present.'));
    
    %% Status bar.
    S.txCurrFile = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [320/W 120/H 175/W 40/H],'string','','backgroundcolor',get(S.fh,'color'),...
        'fontsize',10,'FontUnits','Normalized','HorizontalAlignment','left');
    S.txETA = uicontrol('Parent',S.fh,'style','text','unit','normalized','position',...
        [320/W 80/H 175/W 40/H],'string','','backgroundcolor',get(S.fh,'color'),...
        'fontsize',10,'FontUnits','Normalized','HorizontalAlignment','left');
    
    %% Neuron mode vs. Cardiac mode.
    S.pbMode = uicontrol('Parent',S.fh,'style','push','unit','normalized',...
        'position',[380/W 585/H 110/W 25/H],'FontUnits','Normalized','fontWeight','bold',...
        'string','Cardiac Mode','tooltip',...
        sprintf('Current setting is in cardiac mode.\nClick to change to neuronal mode.'));
    
    %% Run Now! button
    S.pbRunNow = uicontrol('Parent',S.fh,'style','push','unit','normalized',...
        'position',[360/W 20/H 130/W 30/H],'fontsize',11,'FontUnits','Normalized',...
        'string','Start Processing!','callback',{@RunNow,S},'fontWeight','bold',...
        'tooltip',sprintf(['If you have got all the parameters ready,\n' ...
        'click me to Rock-n-Roll !!!']));
    
    %%
    set(S.pbMode,'callback',{@changeMode,S});
    set(S.fh,'ResizeFcn',@figResize);
    %% set gui size gain.
    function figResize(src,evt)     %#ok<INUSD>
        fpos=get(S.fh,'Position');
        FigHeight=round(S.Screen(4)*fpos(4));
        FigWidth=round(FigHeight/H*W);
        fpos(3:4)=[FigWidth/S.Screen(3),FigHeight/S.Screen(4)];
        if fpos(1)<(1-fpos(3))*0.1; fpos(1)=(1-fpos(3))*0.50; end
        set(S.fh,'Position',fpos);
    end
end




function YesNo=isSupportedType(ext)
    supportedType={'.tif';'.tiff';'.lsm';'.lif';'.tf8'};
    YesNo=any(strcmpi(ext,supportedType));
end

%% Add item button.
function [] = AddTiffFile(varargin)
    File=uipickfiles('REFilter', '\.tif$|\.tiff$|\.lif$|\.tf8$|\.lsm$',...
        'Prompt','Select Tiff Image Stacks or Folders');
    
    if isempty(File) || isnumeric(File)
        return
    elseif ischar(File)
        File={File};
    end
    
    S = varargin{3};                    % Get the structure.
    oldstr = get(S.ls,'string');        % The string as it is now.
    
    NewFiles={oldstr{:},File{:}};       %#ok<CCAT>
    set(S.ls,'str',NewFiles,'val',1);
end

%%
function [] = DeleteItem(varargin)
    S = varargin{3};
    List = get(S.ls,'string');
    if ~isempty(List)
        Order = get(S.ls,'value');
        if (Order<1) || (Order>numel(List)); return; end
        List(Order) = [];
        set(S.ls,'string',List,'val',1);
    end
end



function [] = changeMode(varargin)
    S = varargin{3};
    
    str=get(S.pbMode,'string');
    if strcmpi(str,'Cardiac Mode')
        set(S.pbMode,'string','Neuronal Mode','tooltip',...
        sprintf('Current setting is in neuronal mode.\nClick to change to cardiac mode.'));
        set(S.CellMasking,'Value',2);
        set(S.DetectionScaleMin,'String',100);
        set(S.DetectionScaleMax,'String',300);
        set(S.xyt_dim1,'String','0.22');
        set(S.xyt_dim2,'String','0.22');
        set(S.xyt_dim3,'String','50');
        set(S.DetectionLim1,'string','6.75');
        set(S.DetectionLim2,'string','inf');
        set(S.FWHMLim2,'string','100');
        set(S.Tail1,'string','200');
        set(S.Tail2,'string','600');
        set(S.DecayTau,'String','inf');
        set(S.chROIDrawing,'Value',0);
        set(S.chROIApplying,'value',0);
        set(S.chImdilationApplying,'value',1);
        set(S.chWatershed,'Value',0);
    else
        set(S.pbMode,'string','Cardiac Mode','tooltip',...
        sprintf('Current setting is in cardiac mode.\nClick to change to neuronal mode.'));
        set(S.CellMasking,'Value',1);
        set(S.DetectionScaleMin,'String',20);
        set(S.DetectionScaleMax,'String',60);
        set(S.xyt_dim1,'String','0.28');
        set(S.xyt_dim2,'String','0.28');
        set(S.xyt_dim3,'String','10');
        set(S.DetectionLim1,'string','6.27');
        set(S.DetectionLim2,'string','60000');
        set(S.FWHMLim2,'string','10');
        set(S.Tail1,'string','50');
        set(S.Tail2,'string','150');
        set(S.DecayTau,'String','1500');
        set(S.chROIDrawing,'Value',1);
        set(S.chROIApplying,'value',1);
        set(S.chImdilationApplying,'value',0);
        set(S.chWatershed,'Value',1);
    end

end


%% Main code here.
function [] = RunNow(varargin)
    %% Callback for pushbutton, RunNow.
    RunNowButton=varargin{1};
    S = varargin{3};  % Get the structure.
    clear('varargin')

    %% Get all the parameters.
    SparkParameter.xyt_dim = [str2double(get(S.xyt_dim1,'string')),...
        str2double(get(S.xyt_dim2,'string')),...
        str2double(get(S.xyt_dim3,'string'))];
    
    SparkParameter.LocalThreshold=str2num(get(S.LocalThreshold,'string')); %#ok<ST2NM> % NaN for empty
    if isnan(SparkParameter.LocalThreshold);SparkParameter.LocalThreshold=[];end
    if numel(SparkParameter.LocalThreshold)>1;SparkParameter.LocalThreshold=SparkParameter.LocalThreshold(1);end
    
    SparkParameter.GauAmpLim=str2double(get(S.AmpLim,'string'));
    
    SparkParameter.ImageDenoise=(get(S.ImageDenoise,'value'));
    if SparkParameter.ImageDenoise==3
        SparkParameter.ImageDenoise=0;
    end
    
    % % Other parameters.
    DetectionScale=[str2double(get(S.DetectionScaleMin,'string')),NaN,...
                    str2double(get(S.DetectionScaleMax,'string'))];
                
    SparkParameter.DecayTauLim=str2double(get(S.DecayTau,'string'));
    
    SparkParameter.DetectionLim=[str2double(get(S.DetectionLim1,'string')),...
                                 str2double(get(S.DetectionLim2,'string'))];
                             
    SparkParameter.FWHMLim=[str2double(get(S.FWHMLim1,'string')),...
                            str2double(get(S.FWHMLim2,'string'))];
                        
    SparkParameter.SparkSiteMinDistance=str2double(get(S.SparkSiteMinDistance,'string'));
    SparkParameter.DetectorOffset=str2double(get(S.Offset,'string'));
    
    SparkParameter.Leading_Tail=[str2double(get(S.Tail1,'string')),...
                                 str2double(get(S.Tail2,'string'))];
    
    SparkParameter.PhaseShiftCorrect=get(S.Phase,'value')-1;
    SparkParameter.ImRegister=logical(get(S.ImRegister,'value')-1)==1;
    
    SparkParameter.CellMasking=get(S.CellMasking,'value');
    SparkParameter.ROIDrawing=get(S.chROIDrawing,'value')==1;
    SparkParameter.ROIApplying=get(S.chROIApplying,'value')==1;
    SparkParameter.ImdilateApplying=get(S.chImdilationApplying,'value')==1;
    
    SparkParameter.SaveAll=get(S.chSaveAll,'value')==1;
    ImageCorpOption=get(S.chIcrop,'value')==1;
    SparkParameter.UseWatershed=get(S.chWatershed,'Value')==1;

    SparkParameter.RecordPhase=str2double(get(S.RecordPhase,'string'));
    if isnan(SparkParameter.RecordPhase) || SparkParameter.RecordPhase<1
        SparkParameter.RecordPhase=inf;
    end
    
    sCMOSOffset=str2num(get(S.sCMOSdynamicOffsetCorrection,'string')); %#ok<ST2NM>
    sCMOSOffset(2:end)=sCMOSOffset(2:end)/100;
    %% Close user interface.
    set(S.pbAbout,'Enable','off');
    set(S.pbAddFile,'Enable','off');
    set(S.pbDelFolder,'Enable','off');
    % set(RunNowButton,'string','Running!');
    set(S.ls,'Enable','off');
    set(S.xyt_dim1,'Enable','off');
    set(S.xyt_dim2,'Enable','off');
    set(S.xyt_dim3,'Enable','off');
    set(S.LocalThreshold,'Enable','off');
    set(S.AmpLim,'Enable','off');
    set(S.DecayTau,'Enable','off');
    set(S.DetectionLim1,'Enable','off');
    set(S.DetectionLim2,'Enable','off');
    set(S.SparkSiteMinDistance,'Enable','off');
    set(S.FWHMLim1,'Enable','off');
    set(S.FWHMLim2,'Enable','off');
    set(S.sCMOSdynamicOffsetCorrection,'Enable','off');
    set(S.chROIApplying,'Enable','off');
    set(S.chImdilationApplying,'Enable','off');
    set(S.chROIDrawing,'Enable','off');
    set(S.Tail1,'Enable','off');set(S.Tail2,'Enable','off');
    set(S.CellMasking,'Enable','off');
    set(S.chSaveAll,'Enable','off');
    set(S.chShutdown,'Enable','off');
    set(S.ImageDenoise,'Enable','off');
    set(S.DetectionScaleMin,'Enable','off');
    set(S.DetectionScaleMax,'Enable','off');
    set(S.Phase,'Enable','off');
    set(S.ImRegister,'Enable','off');
    set(S.Save,'ButtonDownFcn',{})
    set(S.pbMode,'Enable','off');
    set(S.chIcrop,'Enable','off');
    set(S.Offset,'Enable','off');
    set(RunNowButton, 'Enable', 'off');
    set(RunNowButton,'string','Preparing ...');
    set(S.chWatershed,'Enable','off');
    set(S.RecordPhase,'Enable','off');
    drawnow;
    
    %% Loop all files to get cell mask.
    fprintf('\n\n    +++++++++++++++++++++++  Spark Analysis Started  ++++++++++++++++++++++\n\n\n');
    if (exist('Info_LastDrawingStep2.mat','file')==2)
        LastDrawnFileType=2;
    elseif (exist('Info_LastDrawingStep1.mat','file')==2)
        LastDrawnFileType=1;
    else
        LastDrawnFileType=0;
    end
    
    %  Check whether to load last drawn file.
    if (LastDrawnFileType>0)
        loadLastDrawnInfo=questdlg(['A saved drawn ROI file (Step',num2str(LastDrawnFileType),...
            ') is found. Load it or not?'],'Load Last Drawn?','Yes','No','No');
        loadLastDrawnInfo=strcmpi(loadLastDrawnInfo,'Yes');
    else
        loadLastDrawnInfo=false;
    end
    
    % Check the output file.
    SparkOutFile=get(S.Save,'string');
    if strcmpi(SparkOutFile,'Click to change result saving folder...')            
        SparkOutFileAnswer=questdlg('Output file is not defined. Use the "File-Folder.SparkInfo.txt" as output?',...
            'Spark Info Output Location','Yes','No','Yes');
        SparkOutFileAnswer=strcmpi(SparkOutFileAnswer,'yes');
        if SparkOutFileAnswer
            set(S.Save,'string','File-Folder.SparkInfo.txt');
        else
            return;
        end
    elseif strcmpi(SparkOutFile,'File-Folder.SparkInfo.txt')
        SparkOutFileAnswer=true;
    else
        SparkOutFileAnswer=false;
        [PathName,FileName,Ext]=fileparts(SparkOutFile);
        SparkSiteOutFile=fullfile(PathName,[FileName,'_SparkSites',Ext]);
        clear('PathName','FileName','Ext')
    end
    drawnow expose;

    % Load last drawn or preread all files.
    if loadLastDrawnInfo
        if (LastDrawnFileType==2)
            Info=load('Info_LastDrawingStep2.mat');
        elseif (LastDrawnFileType==1)
            Info=load('Info_LastDrawingStep1.mat');
        end
        Info=Info.Info;
        TotalSeriesNum=numel(Info);
    else        
        % Cycle all the folders and files.
        FolderOrFile = get(S.ls,'string');
        Info=SearchAllFiles(FolderOrFile);  % Info fields: Path,FileName,SeriesNo,TotalNo,Type.
        TotalSeriesNum=numel(Info);

        fprintf('    ========================  Prereading all cells  =======================\n');
        for k=1:TotalSeriesNum
            fprintf('    Step1, reading %d / %d Series, %d%%, %s\n',k,TotalSeriesNum,round(100*k/TotalSeriesNum),Info(k,1).FileName);
            FileName=fullfile(Info(k,1).Path,[Info(k,1).FileName,Info(k,1).Type]);
            [I,I_Info]=ReadStackWithLoci(FileName,Info(k,1).SeriesNo,[1,100]);
            if SparkParameter.PhaseShiftCorrect>0
                I=PhaseShiftAutoAlign(I,SparkParameter.PhaseShiftCorrect);
            end
            
            if any(strcmpi(Info(k,1).Type,{'.tif';'.tiff'}))
                xyt_dim=SparkParameter.xyt_dim;
            else
                xyt_dim=I_Info(1:2,2);
                if xyt_dim(1)<=0; xyt_dim(1)=SparkParameter.xyt_dim(1);end
                if xyt_dim(2)<=0; xyt_dim(2)=SparkParameter.xyt_dim(2);end
            end
            I_Median=mean(I,3);
            % if SparkParameter.CellMasking==4
            %     CurrBW=true(size(I_Median));
            % else
            [CurrBW,~]=CellMasking(I_Median,SparkParameter.CellMasking,xyt_dim);
            % end
            Info(k,1).PreMasking=CurrBW;
            if max(I_Median(:))<=255
                Info(k,1).PreMaskingRaw=uint8(I_Median);
            else
                Info(k,1).PreMaskingRaw=uint16(I_Median);
            end
            clear('I_Median','CurrBW','FileName','I','I_Info')
        end
        save('Info_LastDrawingStep1.mat','Info');
    end
    clear('loadLastDrawnInfo')
    
    %% Manual Checking of cell masks.
    if SparkParameter.ROIDrawing
        if LastDrawnFileType==2
            ROIcheck=questdlg('The pre-drawn ROI Info might have been checked last time, do you really want to check again?',...
                'ROI Checking','Yes','No','Yes');
            ROIcheck=strcmpi(ROIcheck,'yes');
        else
            ROIcheck=true;
        end
        
        if ROIcheck
            fprintf('    ========================  Preparing cell masks  =======================\n');
            for k=1:TotalSeriesNum
                CellMask=FreehandROI(Info(k,1).PreMaskingRaw,Info(k,1).PreMasking,...
                    [num2str(k),' / ',num2str(TotalSeriesNum),' | ',Info(k,1).FileName]);
                CellMask=removePixelROIs(CellMask);
                Info(k,1).PreMasking=CellMask;
                fprintf('    Step2, manual checking %d / %d Series, %d%%, %s\n',k,...
                    TotalSeriesNum,round(100*k/TotalSeriesNum),Info(k,1).FileName);
                clear('CurrBW','CellMask');
            end
            save('Info_LastDrawingStep2.mat','Info');
            fprintf('    ====================  Preparing cell masks finished ===================\n\n\n\n\n\n');
        end
    end
    clear('LastDrawnFileType','ROIcheck')

    %% Loop all files heres.    % Info fields: Path,FileName,SeriesNo,TotalNo,Type.
    ElapsedTime=tic;
    set(RunNowButton,'string','Running ...');
    
    for k=1:TotalSeriesNum
        % Print current file information in the command window.
        fprintf('    Overall progress :    %d / %d Series, %d%%\n\n',...
            k,TotalSeriesNum,round(100*k/TotalSeriesNum));
        fprintf('    Current Location :    %s\n',OrganizePath(Info(k,1).Path,46));
        fprintf('    Current File Name:    %s\n',[Info(k,1).FileName Info(k,1).Type]);
        if Info(k,1).TotalNo>1
            fprintf('    Current Sereis   :    %0.0f / %0.0f\n',Info(k,1).SeriesNo,Info(k,1).TotalNo);
        end
        % Progress bar.
        set(S.txCurrFile,'string',sprintf(['Current Series: No. ',num2str(k),' / ',...
            num2str(TotalSeriesNum),', ',num2str(round(k/TotalSeriesNum*100)),'%%']));
        drawnow expose;
        
        FileName=fullfile(Info(k,1).Path,[Info(k,1).FileName,Info(k,1).Type]);
        
        % %%%%%%%%%%%%%%%%%%%% Real Analysis Start Here %%%%%%%%%%%%%%%%%%%%%%%%
        try                    % Use try-catch control to escape from error.
            % Read image series.
            fprintf('    Reading file from Hard Disk ...    ');
            [I,I_Info]=ReadStackWithLoci(FileName,Info(k,1).SeriesNo);
            if I_Info(5,1)>1; I=I(:,:,1:I_Info(5,1):end); end % Read the first channel to analyze.
            
            Isiz=size(I);
            
            if Isiz(3)<40    % Check whether the image series is just too small.
                fprintf('\n    This file does not look like a spark movie, skip it now.\n\n');
                RemainedTime=toc(ElapsedTime)/k*(TotalSeriesNum-k);
                set(S.txETA,'string',sprintf(['Remaining Time:\n',num2str(floor(RemainedTime/3600)),...
                    ' h ',num2str(ceil(mod(RemainedTime,3600)/60)),' min']));
                drawnow expose;
                clear('I','FileName','SparkLogFileName','SparkSiteLogName','I_Info',...
                    'StackFileName','RemainedTime');
                continue;
            end
            fprintf('    done\n');
            
            % Correct sCMOS dynamic
            if ~(isnan(sCMOSOffset(1)) || (round(sCMOSOffset(1))>2) || (round(sCMOSOffset(1))<1))
                fprintf('    sCMOS Dynamic Offset correction ...');
                I=sCMOSdynamicOffsetCorrection(I,sCMOSOffset);
                fprintf('    done\n');
            end
            
            % Crop unnecessary regions to accelerate the processing.
            CellROI=bwlabeln(Info(k,1).PreMasking);
            numCellROI=max(CellROI(:));
            for ROIno=1:numCellROI
                CellMask=CellROI==ROIno;
                if ImageCorpOption
                    fprintf('    Cropping image stack ...       ');
                    bwCrop=regionprops(CellMask,'BoundingBox');
                    ROInum=numel(bwCrop);
                    x1=zeros(ROInum,1); x2=zeros(ROInum,1); y1=zeros(ROInum,1); y2=zeros(ROInum,1);
                    for j=1:numel(bwCrop)
                        x1(j)=floor(bwCrop(j).BoundingBox(2));    if x1(j)<1; x1(j)=1; end
                        x2(j)=x1(j)+bwCrop(j).BoundingBox(4)+1;   if x2(j)>Isiz(1); x2(j)=Isiz(1); end
                        y1(j)=floor(bwCrop(j).BoundingBox(1));    if y1(j)<1; y1(j)=1; end
                        y2(j)=y1(j)+bwCrop(j).BoundingBox(3)+1;   if y2(j)>Isiz(2); y2(j)=Isiz(2); end
                    end
                    x1=min(x1); x2=max(x2); y1=min(y1); y2=max(y2);
                    
                    I=I(x1:x2,y1:y2,:);
                   CellMask=CellMask(x1:x2,y1:y2);
                    clear('bwCrop','ROInum','j');
                    fprintf('        done\n');
                end
                
                
                if any(strcmpi(Info(k,1).Type,{'.tif';'.tiff'}))
                    xyt_dim=SparkParameter.xyt_dim;
                else
                    xyt_dim=I_Info(1:3,2); xyt_dim(3)=xyt_dim(3)*1000;
                    if xyt_dim(1)<=0; xyt_dim(1)=SparkParameter.xyt_dim(1);end
                    if xyt_dim(2)<=0; xyt_dim(2)=SparkParameter.xyt_dim(2);end
                    if xyt_dim(3)<=0; xyt_dim(3)=SparkParameter.xyt_dim(3);end
                end
                
                
                % Real analysis of the image series.
                DetectionScale(2)=max(1,round((DetectionScale(3)-DetectionScale(1))/5/xyt_dim(3)))*xyt_dim(3);
                if isinf(SparkParameter.RecordPhase) || (SparkParameter.RecordPhase>=Isiz(3))
                    if numCellROI>1
                        fprintf('\n    Processing Recording ROI %d:\n',ROIno);
                    end
                    SparkAll=SparkAnalysis(I,...
                        'xyt_dim',              xyt_dim,...
                        'Threshold',            SparkParameter.LocalThreshold,...
                        'CameraOffset',         SparkParameter.DetectorOffset,...
                        'GauAmpLim',            SparkParameter.GauAmpLim,...
                        'TauLim',               SparkParameter.DecayTauLim,...
                        'DetectionLimit',       SparkParameter.DetectionLim,...
                        'FWHMLimit',            SparkParameter.FWHMLim,...
                        'SparkLim',             SparkParameter.Leading_Tail,...
                        'Denoising',            SparkParameter.ImageDenoise,...
                        'DetectionScale',       DetectionScale,...
                        'MinSparkSiteDistance', SparkParameter.SparkSiteMinDistance,...
                        'CellMasking',          CellMask,...
                        'PhaseShiftCorrect',    SparkParameter.PhaseShiftCorrect,...
                        'ImRegister',           SparkParameter.ImRegister,...
                        'ApplyGlobalCellmask',  SparkParameter.ROIApplying,...
                        'ApplyImdilation',      SparkParameter.ImdilateApplying,...
                        'UseWatershedSegment',  SparkParameter.UseWatershed);
                else
                    numRecord=floor(Isiz(3)/SparkParameter.RecordPhase);
                    % numRecord=2; % for testing
                    SparkAll=[];
                    for recordphaseNo=1:numRecord
                        if numCellROI==1
                            fprintf('\n    Processing Recording Phase No.: %d / %d\n',recordphaseNo,numRecord);
                        else
                            fprintf('\n    Processing ROI %d, Recording Phase No.: %d / %d\n',ROIno,recordphaseNo,numRecord);
                        end
                        z1=1+(recordphaseNo-1)*SparkParameter.RecordPhase;
                        z2=recordphaseNo*SparkParameter.RecordPhase;
                        SparkAllTemp=SparkAnalysis(I(:,:,z1:z2),...
                            'xyt_dim',              xyt_dim,...
                            'Threshold',            SparkParameter.LocalThreshold,...
                            'CameraOffset',         SparkParameter.DetectorOffset,...
                            'GauAmpLim',            SparkParameter.GauAmpLim,...
                            'TauLim',               SparkParameter.DecayTauLim,...
                            'DetectionLimit',       SparkParameter.DetectionLim,...
                            'FWHMLimit',            SparkParameter.FWHMLim,...
                            'SparkLim',             SparkParameter.Leading_Tail,...
                            'Denoising',            SparkParameter.ImageDenoise,...
                            'DetectionScale',       DetectionScale,...
                            'MinSparkSiteDistance', SparkParameter.SparkSiteMinDistance,...
                            'CellMasking',          CellMask,...
                            'PhaseShiftCorrect',    SparkParameter.PhaseShiftCorrect,...
                            'ImRegister',           SparkParameter.ImRegister,...
                            'ApplyGlobalCellmask',  SparkParameter.ROIApplying);
                        SparkAll=cat(1,SparkAll,SparkAllTemp);
                        clear('SparkAllTemp','z1','z2')
                    end
                    clear('numRecord','recordphaseNo');
                end

                % Write metainformation to the result file
                if Info(k,1).TotalNo==1
                    StackFileName=FileName;
                elseif Info(k,1).TotalNo>1
                    StackFileName=[FileName,'_Series_',num2str(Info(k,1).SeriesNo)];
                end
                
                if numCellROI>1
                    StackFileName=[StackFileName,'_ROI',num2str(ROIno)]; %#ok<AGROW>
                end
                
                if SparkOutFileAnswer
                    SparkOutFile=fullfile(Info(k,1).Path,'SparkInfo.txt');
                    SparkSiteOutFile=fullfile(Info(k,1).Path,'SparkSitesInfo.txt');
                end
                
                fprintf('    Saving Spark/Site results to Hard Disk...      ');
                saveTxtResults(SparkAll,StackFileName,SparkOutFile,SparkSiteOutFile);
                fprintf('  done\n');
                
                % Save all intermediate variables.
                if SparkParameter.SaveAll
                    fprintf('    Saving all intermediate results to Hard Disk...');
                    % SparkAll.FileName=StackFileName;
                    [~,SparkAll(1).FileName,~]=fileparts([StackFileName,'.mat']);
                    save([StackFileName,'.mat'],'SparkAll','-v7.3');
                    fprintf('  done\n')
                end
            end
            % Progress bar.
            RemainedTime=toc(ElapsedTime)/k*(TotalSeriesNum-k);
            set(S.txETA,'string',sprintf(['Remaining Time:\n',num2str(floor(RemainedTime/3600)),...
                ' h ',num2str(ceil(mod(RemainedTime,3600)/60)),' min']));
            drawnow expose;
            
            % Clear the resulted variables.
            clear('SparkAll','I','FileName','SparkLogFileName','SparkSiteLogName','I_Info',...
                'xyt_dim','StackFileName','RemainedTime');
        catch err
            disp(err.message);
        end
        % %%%%%%%%%%%%%%%%%%%% Analysis End Here     %%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n\n\n\n\n\n');
    end
    
    fprintf('    +++++++++++++++++++++++ Spark Analysis Finihsed ++++++++++++++++++++++\n\n\n\n');
    
    %%  Shutdown the system, if possible.
    if get(S.chShutdown,'value')==1
        try
            if ispc;    eval('!shutdown -s -f -t 60');  quit;   end
            if isunix;  eval('!shutdown -h +1');  quit; end
        catch err
            disp(err.message);
        end
    end
    %% Set gui available again.
    set(S.pbAbout,'Enable','on');
    set(S.pbAddFile,'Enable','on');
    set(S.pbDelFolder,'Enable','on');
    set(RunNowButton, 'Enable', 'on');
    set(RunNowButton,'string','Start Processing!');
    set(S.ls,'Enable','on');
    set(S.xyt_dim1,'Enable','on');
    set(S.xyt_dim2,'Enable','on');
    set(S.xyt_dim3,'Enable','on');
    set(S.ImRegister,'Enable','on');
    set(S.LocalThreshold,'Enable','on');
    set(S.AmpLim,'Enable','on');
    set(S.DecayTau,'Enable','on');
    set(S.DetectionLim1,'Enable','on');
    set(S.DetectionLim2,'Enable','on');
    set(S.SparkSiteMinDistance,'Enable','on');
    set(S.FWHMLim1,'Enable','on');
    set(S.FWHMLim2,'Enable','on');
    set(S.sCMOSdynamicOffsetCorrection,'Enable','on');
    set(S.chROIApplying,'Enable','on');
    set(S.chImdilationApplying,'Enable','on');
    set(S.chROIDrawing,'Enable','on');
    set(S.Tail1,'Enable','on');set(S.Tail2,'Enable','on');
    set(S.CellMasking,'Enable','on');
    set(S.chSaveAll,'Enable','on');
    set(S.ImageDenoise,'Enable','on');
    set(S.DetectionScaleMin,'Enable','on');
    set(S.DetectionScaleMax,'Enable','on');
    set(S.chShutdown,'Enable','on');
    set(S.Offset,'Enable','on');
    set(S.txCurrFile,'string','');
    set(S.txETA,'string','');
    set(S.Save,'ButtonDownFcn',{@CreateSaveFile,S})
    set(S.Phase,'Enable','on');
    set(S.pbMode,'Enable','on');
    set(S.chIcrop,'Enable','on');
    set(S.chWatershed,'Enable','on');
    set(S.RecordPhase,'Enable','on');
    drawnow expose;
end

function CellMask=removePixelROIs(CellMask)
    CellMask=imclose(CellMask,strel('disk',1));
    CellMask=imopen(CellMask,strel('disk',1));
    CellMask=imfill(CellMask,'holes');
    
    % Find out the largest area to serve as a cell.
    CellMaskLabel = bwconncomp(CellMask,ones(3,3));
    area = cellfun(@numel, CellMaskLabel.PixelIdxList);
    idxToKeep=area>=max(area*0.05);
    CellMaskLabel.PixelIdxList=CellMaskLabel.PixelIdxList(idxToKeep);
    CellMaskLabel.NumObjects=sum(idxToKeep);
    
    CellMask=false(size(CellMask));
    for k = 1 : CellMaskLabel.NumObjects
        CellMask(CellMaskLabel.PixelIdxList{k}) = true;
    end
end

%% Search all files supported.
function Info=SearchAllFiles(FolderOrFile)
    Info=struct('Path',{},'FileName',{},'SeriesNo',{},'TotalNo',{},'Type',{});
    for k=1:numel(FolderOrFile)
        if exist(FolderOrFile{k},'dir')==7
            CurrList=dir(fullfile(FolderOrFile{k}));
            for j=1:numel(CurrList)
                if ~strcmpi(CurrList(j).name(1),'.')
                    Info=cat(1,Info,SearchAllFiles({fullfile(FolderOrFile{k},CurrList(j).name)}));
                end
            end
        elseif exist(FolderOrFile{k},'file')==2
            [CurrPath,CurrFileName,CurrExt]=fileparts(FolderOrFile{k});
            if isSupportedType(CurrExt)
                try
                    if strcmpi(CurrExt,'.tif') || strcmpi(CurrExt,'.tiff')
                        numSeries=1;
                    else
                        numSeries=ReadStackWithLoci(FolderOrFile{k},'GetStackNum');
                    end
                    m=numel(Info);
                    for j=1:numSeries
                        Info(m+j,1).Path=CurrPath;
                        Info(m+j,1).FileName=CurrFileName;
                        Info(m+j,1).SeriesNo=j;
                        Info(m+j,1).TotalNo=numSeries;
                        Info(m+j,1).Type=CurrExt;
                    end
                catch err
                    disp(err.message);
                end
            end
        end
    end
end


function [] = CreateSaveFile(varargin)
    S = varargin{3};                    % Get the structure.
    [FileName,PathName]=uiputfile('.txt','Save spark info to','SparkInfo');
    if ~isscalar(FileName)
        set(S.Save,'string',fullfile(PathName,FileName));
    else
        return
    end
end
%%
function str=OrganizePath(str,strlenToKeep)
    str_len=length(str);
    if str_len>strlenToKeep; str=['...',str((str_len-(strlenToKeep-1)):end)]; end
end


%%
function [] = About(varargin)
    msgbox({'iSpark (interlligent Spark analysis) can do:';...
        '    1. PURE-LET to remove Poissonian noise.';...
        '    2. Automatic spark detection.';...
        '    3. Automatic spark validation.';...
        '    4. Batch analysis of calcium spark properties.';...
        '';...
        'Precedures to do a typical analysis:';...
        '    1. Click "Add " to add single image stack or stack collection files to analyze;';...
        '    2. Set the pixel dimensions (x/y/t) in um and ms;';...
        '    3. Click "Start Processing!" to start the analysis;';...
        '';...
        ['Advanced parameters are also changable.' ...
        'Type "help SparkAnalysis" for more information.'];...
        '';...
        'Quick help: move the mouse over a parameter to get quick help info.';...
        '';...
        ['Note: all the analysis results will be automatically saved into text ', ...
        'files specified in the output dialog, including the information of', ...
        'single sparks and spark sites.'];...
        '';...
        'Developed by Qinghai Tian (www.lipplab.de).';...
        '';...
        'Affiliation:';...
        '    Institute for Molecular Cell Biology, Research Center for Molecular';...
        '    Imaging and Screening. Medical Faculty, Saarland University.';...
        '    Building 61, 66421 Homburg/Saar, Germany.'},'About','help');
end

%%
function [Pos,H,W,Screen]=FigSize
    H=680;W=510;    Screen=get(0,'ScreenSize');
    FigHeight=round(Screen(4)*0.668);    if FigHeight<H; FigHeight=H; end
    FigWidth=round(FigHeight/H*W);
    Pos=[round((Screen(3)-FigWidth)*0.5),round((Screen(4)-FigHeight)*0.7),FigWidth,FigHeight];
    Pos(1)=Pos(1)/Screen(3); Pos(3)=Pos(3)/Screen(3);
    Pos(2)=Pos(2)/Screen(4); Pos(4)=Pos(4)/Screen(4);
end
