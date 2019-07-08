function iSparkControl(SparkAll,IdToShow)
    %iSparkControl, a GUI-based single spark contol utility.
    %
    %Input:
    %    SparkAll, the struct that is returned from the main spark program.
    %
    %Usage:
    %    iSparkControl;
    %          or
    %    iSparkControl(SparkAll);
    %          or
    %    iSparkControl(SparkAll,ID); % just to show single spark: 
    
    
    %% Plot all the sparks indicated in the second input.
    if nargin==2
        if isfield(SparkAll,'Creator')
            if ~isempty(strfind(SparkAll.Creator,'park'))
                if ~isfloat(SparkAll.Data)
                    SparkAll.Data=double(SparkAll.Data-SparkAll.CameraOffset);
                end
                SparkAll.Bgr=mean(SparkAll.Data,3);
                IDmax=size(SparkAll.SparkPos,1);
            end
        end
        
        if isscalar(IdToShow)
            IdToShow=round(IdToShow);
            for k=1:numel(IdToShow)
                if (IdToShow(k)>0) && (IdToShow(k)<=IDmax)
                    Live_fig=figure;
                    set(Live_fig,'Position',SparkFigSize);
                    figure(Live_fig);
                    ShowSpark(SparkAll,Live_fig,IdToShow(k),false,[]);
                end
            end
            return;
        elseif ischar(IdToShow) && strcmpi(IdToShow,'rejected')
            SparkAll.SpatiotemporalFitting=SparkAll.SpatiotemporalFittingRejected;
            SparkAll.SparkProperty=SparkAll.SparkPropertyRejected;
            SparkAll.SparkPos=SparkAll.SparkPosRejected;
        end
        
    end
    
    %%
    GuiSize=FigSize;
    S.fh = figure('units','pixels','position',GuiSize,'menubar','none','name','iSpark Control',...
        'numbertitle','off','resize','off','Color',[0.95 0.95 0.95],'HandleVisibility','off');
    %% Input
    if nargin>=1
        if isfield(SparkAll,'Creator')
            if ~isempty(strfind(SparkAll.Creator,'park'))
                if ~isfloat(SparkAll.Data)
                    SparkAll.Data=single(SparkAll.Data)-SparkAll.CameraOffset;
                end
                SparkAll.Bgr=mean(SparkAll.Data,3);
                setappdata(S.fh,'SparkAll',SparkAll);
                setappdata(S.fh,'File','from workspace');
                setappdata(S.fh,'Path','from workspace');
                setappdata(S.fh,'LiveFigureHandle',[]);
                IDmax=size(SparkAll.SparkPos,1);
            else
                IDmax=0;
            end
        else
            IDmax=0;
        end
        
        if nargin<2; clear('SparkAll');end
        if IDmax==0
            msgbox('The input is not correct. Please load it again.','Incorrect Input','warn');
        end
    else
        IDmax=0;
    end
    
    %% Load file.
    S.pbLoadFile = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[5 GuiSize(4)-35 57 30],'string','Load ...',...
        'fontWeight','bold','tooltip','Click to load saved iSPark intermediate results');
    %% Disp S.
    S.pbDisplayContent = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[65 GuiSize(4)-35 57 30],'string','Metadata',...
        'fontWeight','bold','tooltip','Display metadata of current results');
    
    %% Next Buttons.
    S.pbPrevious = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[140 GuiSize(4)-35 25 30],'string','<','fontsize',14,...
        'fontWeight','bold','tooltip','Display previous spark');
    S.ID = uicontrol('Parent',S.fh,'style','edit','unit','pix','position',...
        [168 GuiSize(4)-33 50 27],'string','','backgroundcolor',get(S.fh,'color'),...
        'tooltip','Spark IDs to display.','fontWeight','bold','fontsize',12);
    S.pbNext = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[221 GuiSize(4)-35 25 30],'string','>','fontsize',14,...
        'fontWeight','bold','tooltip','Display next spark');
    %% Manual check.
    S.pbManualCheck = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[264 GuiSize(4)-35 85 30],'string','Manaul Check',...
        'fontWeight','bold','tooltip',sprintf(['Check the sparks indicated in\n' ...
        'the Spark ID box one by one...\n \n' ...
        'To view all the sparks, leave the\nSpark ID box EMPTY.']));
    %% LiveView.
    S.pbLiveView = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[352 GuiSize(4)-35 70 30],'string','Live View',...
        'fontWeight','bold','tooltip',sprintf(['Visualize the sparks indicated in' ...
        ' the Spark\nID box, automatically and one by one...\n \n' ...
        'To view all the sparks, leave the\nSpark ID box EMPTY.']));
    %% Set special tooltips.
    if IDmax>=1
        set(S.ID,'tooltip',sprintf(['Spark IDs to display.\n \nThe range is between 1 ~ ' ...
            num2str(IDmax) '.']));
    else
        set(S.ID,'tooltip',sprintf(['Spark IDs to display.\n \n' ...
            'Currently no result is loaded.\nPlease load it now.']));
    end
    
    %% Set callback.
    set(S.pbLoadFile,       'callback',{@LoadSavedMatFile,S});
    set(S.pbDisplayContent, 'callback',{@DispS,S});
    
    set(S.pbLiveView,   'callback',{@LiveView,S});
    set(S.pbManualCheck,'callback',{@LiveView,S});
    
    set(S.ID,        'callback',{@PreNextSpark,S});
    set(S.pbPrevious,'callback',{@PreNextSpark,S});
    set(S.pbNext,    'callback',{@PreNextSpark,S});
    
end













%% Gui Position.
function Pos=FigSize
    Screen=get(0,'ScreenSize');
    FigHeight=39;
    FigWidth=425;
    Pos=[10 round((Screen(4)-FigHeight-30)) FigWidth FigHeight];
end

%% Add item button.
function LoadSavedMatFile(varargin)
    S = varargin{3};                    % Get the structure.
    Path=getappdata(S.fh,'Path');
    [File,Path] = uigetfile(fullfile(Path,'*.mat'),'Select a saved iSpark intermediate result file');
    if File==0; return; end
    R=load(fullfile(Path,File));
    SparkAll=fieldnames(R);
    set(S.ID,'tooltip',sprintf(['Spark IDs to display.\n \nThe range is between 1 ~ ' ...
        num2str(size(R.(SparkAll{1}).SparkPos,1)) '.']));
    if ~isfloat(R.(SparkAll{1}).Data)
        R.(SparkAll{1}).Data=double(R.(SparkAll{1}).Data)-R.(SparkAll{1}).CameraOffset;
    end
    R.(SparkAll{1}).Bgr=mean(R.(SparkAll{1}).Data,3);
    setappdata(S.fh,'SparkAll',R.(SparkAll{1}));
    setappdata(S.fh,'File',File);
    setappdata(S.fh,'Path',Path);
end

%% Disp struct.
function DispS(varargin)
    S_fig = varargin{3};                    % Get the structure.
    S=getappdata(S_fig.fh,'SparkAll');
    if isempty(S)
        msgbox('No result is loaded. Please load it now.','No Result','warn');
        return
    end
    Path=getappdata(S_fig.fh,'Path');
    if strcmpi(Path(end),filesep); Path=Path(1:(numel(Path)-1)); end
    [~,File,Ext]=fileparts(Path);    Path=['...',filesep,File,Ext,filesep];
    str1=sprintf('    Currrent Path:    %s\n    Currrent File:    %s',...
        Path,getappdata(S_fig.fh,'File'));
    if isfield(S,'SparkPos')
        SparkNum=size(S.SparkPos,1);
    else
        SparkNum=0;
    end
    
    if isfield(S,'SparkSitePos')
        SparkSiteNum=size(S.SparkSitePos,1);
    else
        SparkSiteNum=0;
    end
    
    strsp=newline;
    str2=sprintf('    Spark Number:       %d',SparkNum);
    str3=sprintf('    Spark Site Number:  %d',SparkSiteNum);
    str4=sprintf('    Device Gain:   %0.2f',S.Gain(1));
    str5=sprintf('    Smoothn s:     %0.4f',S.Gain(2));
    str6=sprintf('    Pixel Dimension:  %0.2fum x %0.2fum x %0.1fms',...
        S.xyt_dim(1),S.xyt_dim(2),S.xyt_dim(3));
    str7=sprintf('    File Size:     %d x %d x %d pixels',...
        size(S.Data,1),size(S.Data,2),size(S.Data,3));
    str8=sprintf('    Cell Area:     %0.2f sqare um',...
        sum(double(S.CellMask(:)))*S.xyt_dim(1)*S.xyt_dim(2));
    strtitle=sprintf('Current Metadata:');
    msgbox({strtitle;strsp;str1;strsp;str7;str6;str8;strsp;...
        str4;str5;strsp;str2;str3;strsp},'Current Metadata','help');
end



%% LiveView
function LiveView(varargin)
    S = varargin{3};                    % Get the structure.
    if varargin{1}==S.pbManualCheck
        LiveShow=false;
    elseif varargin{1}==S.pbLiveView
        LiveShow=true;
    else
        LiveShow=true;
    end
    
    
    % Get the contents to show.
    if ~isappdata(S.fh,'SparkAll')
        msgbox('No result is loaded. Please load it now.','No Result','warn');
        return
    else
        SparkAll=getappdata(S.fh,'SparkAll');
    end
    
    % Get the ID to show.
    ID=str2num(get(S.ID,'string'));     %#ok<*ST2NM>
    if isempty(ID)
        ID=1:size(SparkAll.SparkPos,1);
    else
        ID=ID((ID>=1) & (ID<=size(SparkAll.SparkPos,1)));
    end
    if numel(ID)==0
        msgbox('No valid spark ID. Please update it.','No Result','warn');
        return
    end
    
    % Get pause time for every spark.
    if LiveShow
        PauseTime=inputdlg(sprintf(['Live View:    <   %d   >   sparks to show.' ...
            '\n \nInput a pause time for every spark in second:\n'],numel(ID)),...
            'Pause Time',1,{'0.3'});
        if isempty(PauseTime)
            return;
        else
            PauseTime=abs(str2num(PauseTime{1}));
        end
    end
    
    % Get the live figure handle.
    Live_fig=getappdata(S.fh,'LiveFigureHandle');
    if isempty(Live_fig)
        Live_fig=figure('resize','off','menubar','none','toolbar','figure',...
            'name','Live Figure','numbertitle','off');
        set(Live_fig,'Position',SparkFigSize); figure(S.fh);
        setappdata(S.fh,'LiveFigureHandle',Live_fig);
        KillFig=true;
    elseif ishandle(Live_fig)==0
        Live_fig=figure('resize','off','menubar','none','toolbar','figure',...
            'name','Live Figure','numbertitle','off');
        set(Live_fig,'Position',SparkFigSize); figure(S.fh); % set(Live_fig,'PaperUnits','points');
        setappdata(S.fh,'LiveFigureHandle',Live_fig);
        KillFig=true;
    else
        KillFig=false;
    end
    
    
    % Loop all sparks.
    if ~LiveShow
        msgtitle=['    No.\tVar(cx,cy,ct)\tMass(squm*ms)\tFWHM(um)\tAxisMax(um)\tAxisMin(um)'... 
        '\tAmpFFzero\tTau(ms)\tUpstrokeSigma(ms)\tAmpInt\tBgrInt\tDistToMem(um)\n'];
        fprintf(msgtitle);
    end
    
    LastPos=[];
    Accept=0;
    Reject=0;
    for k=1:numel(ID)
        [StatusGood,LastPos]=ShowSpark(SparkAll,Live_fig,ID(k),~LiveShow,LastPos);
        if strcmp(StatusGood,'Accept')
            Accept=Accept+1;
        elseif strcmp(StatusGood,'Reject')
            Reject=Reject+1;
        else
            break;
        end
        if LiveShow
            drawnow;
            pause(PauseTime);
        end
        % % Save as pdf.
        % set(Live_fig,'PaperOrientation','landscape');
        % set(Live_fig,'position',[100,100,900,450]);
        % print(Live_fig,'iSparkDetectionExample.ps','-dpsc','-append');
        % print(Live_fig,['iSparkDetectionExample_',num2str(ID(k)),'.ps'],'-dpsc','-append')
        % saveas(Live_fig,['Example_',num2str(ID(k)),'.pdf']);
    end
    if ~LiveShow
        fprintf('\n\n\tAccepted: %g;  Rejected: %g;  Total: %g\n',Accept,Reject,Accept+Reject);
    end
    
    % Close the live figure.
    if KillFig && ishandle(Live_fig)
        close(Live_fig);
    end
end




%% PreviousSpark button
function PreNextSpark(varargin)
    
    S = varargin{3};                    % Get the structure.
    if varargin{1}==S.pbPrevious
        PreNext=true;
        SStep=-1;
    elseif varargin{1}==S.pbNext
        PreNext=true;
        SStep=1;
    elseif varargin{1}==S.ID
        PreNext=false;
        SStep=0;
    end
    
    % Get the contents to show.
    if ~isappdata(S.fh,'SparkAll')
        msgbox('No result is loaded. Please load it now.','No Result','warn');
        return
    end
    
    % Get the ID to show.
    ID=str2num(get(S.ID,'string'));     %#ok<*ST2NM>
    if isempty(ID)
        if PreNext
            msgbox('Spark ID is empty. Please input it now.','Spark ID','warn');
        end
        return
    elseif numel(ID)>1
        if PreNext
            msgbox('More than one ID found. Please input single ID now.','Spark ID','warn');
        end
        return
    end
    
    SparkAll=getappdata(S.fh,'SparkAll');
    % Last spark
    if (SStep==1) && (ID==size(SparkAll.SparkPos,1))
        msgbox('Current spark is already the last spark.','Spark ID','warn');
        return
    end
    
    % first spark
    if (SStep==-1) && (ID==1)
        msgbox('Current spark is already the first spark.','Spark ID','warn');
        return
    end
    
    % Out of range
    if(ID<1) || (ID>size(SparkAll.SparkPos,1))
        msgbox('Spark ID is out of range. Please correct it now.','Spark ID','warn');
        return
    end
    
    ID=ID+SStep;
    set(S.ID,'string',num2str(ID));
    
    % Get the live figure handle.
    Live_fig=getappdata(S.fh,'LiveFigureHandle');
    if isempty(Live_fig)
        Live_fig=figure('resize','off','menubar','none','toolbar','figure',...
            'name','Live Figure','numbertitle','off');
        set(Live_fig,'Position',SparkFigSize); figure(S.fh);
        setappdata(S.fh,'LiveFigureHandle',Live_fig);
    elseif ishandle(Live_fig)==0
        Live_fig=figure('resize','off','menubar','none','toolbar','figure',...
            'name','Live Figure','numbertitle','off');
        set(Live_fig,'Position',SparkFigSize); figure(S.fh);
        setappdata(S.fh,'LiveFigureHandle',Live_fig);
    end
    % set(Live_fig,'PaperUnits','points');
    
    ShowSpark(SparkAll,Live_fig,ID,false,[]);
    
end



%%
function [StatusGood,MFquestdlgPos]=ShowSpark(I,Live_fig,ID,ManualCheck,MFquestdlgPos)
    StatusGood='Reject';
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    end
    set(Live_fig,'Name',['Single Spark No. ' num2str(ID)],'NumberTitle','Off')
    RawSize=size(I.Data);

    x=round(I.SparkPos(ID,7));
    y=round(I.SparkPos(ID,8));
    t=I.SparkPos(ID,9);
    if x<1;x=1;end; if x>RawSize(1);x=RawSize(1);end
    if y<1;y=1;end; if y>RawSize(2);y=RawSize(2);end
    if t<1;t=1;end; if t>RawSize(3);t=RawSize(3);end
    
    tracefull=squeeze(I.Data(x,y,:));
    trace=I.SpatiotemporalFitting{ID,2};
    
    x1=I.SparkPos(ID,1);x2=I.SparkPos(ID,2);
    y1=I.SparkPos(ID,3);y2=I.SparkPos(ID,4);
    t1=I.SparkPos(ID,5);t2=I.SparkPos(ID,6);

    % Time projection the spark.
    sparkTimeProject=I.SpatiotemporalFitting{ID,3};
    
    % Spark montage.
    sparkMontage=I.Data(x1:x2,y1:y2,t1:t2);
    sparkMontage=bsxfun(@minus,sparkMontage,I.Bgr(x1:x2,y1:y2));
    
    tFrame=sum(sum(sparkMontage,1),2);
    tFrame=tFrame(:);
    tFrame=find(tFrame==max(tFrame));
    tFrame=tFrame(1);
    sparkGlobalPos=I.Data(:,:,ceil(t1+tFrame-1));
    sparkGlobalPos(x1:x2,y1)=max(sparkGlobalPos(:));
    sparkGlobalPos(x1:x2,y2)=max(sparkGlobalPos(:));
    sparkGlobalPos(x1,y1:y2)=max(sparkGlobalPos(:));
    sparkGlobalPos(x2,y1:y2)=max(sparkGlobalPos(:));
    
    

    LowLim=sort(sparkMontage(:),'ascend');
    if numel(LowLim)<4 || numel(LowLim)>0
        LowLim=LowLim(1);
    elseif numel(LowLim)>=4
        LowLim=LowLim(round(numel(LowLim)*0.2));
    else
        msgbox(['Spark ' num2str(ID) 'is not a real spak.'],'No Result','error');
    end

    sparkMontage=sparkMontage(:,:,1:max(1,round(6/I.xyt_dim(3))):end);
    sparkMontage=reshape(sparkMontage,[x2-x1+1,y2-y1+1,1,size(sparkMontage,3)]);
    [sparkMontage,space_bw]=montage_Modified(sparkMontage,[NaN floor(sqrt(size(sparkMontage,4)))]);
    sparkMontage(space_bw)=LowLim;
    
    clear('x','y','tt2','x1','x2','y1','y2','xx1','xx2','yy1','yy2','space_bw','temporalExtend')
    
    sparkAmp=I.SparkProperty(ID,4);
    sparkBgr=I.SparkProperty(ID,5);
    sparkFF0=sparkAmp/sparkBgr;
    
    %% Global view.
    if size(sparkGlobalPos,1)>size(sparkGlobalPos,2);sparkGlobalPos=sparkGlobalPos';end
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        GlobalViewH=subplot(2,20,1:12,'Parent',Live_fig);
    end
    imagesc(sparkGlobalPos,'Parent',GlobalViewH,[min(sparkGlobalPos(:)) max(sparkGlobalPos(:))]);
    axis(GlobalViewH,'image','off'); 
    if size(sparkGlobalPos,1)/size(sparkGlobalPos,2)>0.5
        colorbar('peer',GlobalViewH,'EastOutside');
    else
        colorbar('peer',GlobalViewH,'SouthOutside');
    end

    title(GlobalViewH,['Spark No.' num2str(ID) '    |    dAmp=' ...
        num2str(round(sparkAmp)) '    |    Bgr=' num2str(round(sparkBgr)) ...
        '    |    dF/F0=' num2str(round(sparkFF0*100)/100)]);
    colormap(GlobalViewH,coolcolor);
    
    %% 2D view.
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        View2DH=subplot(2,20,[17 20],'Parent',Live_fig);
    end
    try
        imagesc(sparkTimeProject,'Parent',View2DH,[LowLim max(sparkTimeProject(:))]);
    catch
        imagesc(sparkTimeProject,'Parent',View2DH);
    end
    axis(View2DH,'image','off');
    title(View2DH,'Temporal projection');
    colormap(View2DH,coolcolor);
    
    %% 2D Montage View.
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        MontageH=subplot(2,20,[13 16],'Parent',Live_fig);
    end
    try
        imagesc(sparkMontage,'Parent',MontageH,[LowLim,max(sparkMontage(:))]);
    catch
        imagesc(sparkMontage,'Parent',MontageH);
    end
    colorbar('peer',MontageH,'EastOutside');
    title(MontageH,'Spark Montage');
    axis(MontageH,'off','image');
    colormap(MontageH,coolcolor);
    %% Single Spark trace.
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        SingleTraceH=subplot(2,20,[37 40],'Parent',Live_fig);
    end
    plot(SingleTraceH,trace(:,1)*I.xyt_dim(3),trace(:,2),'o');
    hold(SingleTraceH,'all');
    plot(SingleTraceH,trace(:,1)*I.xyt_dim(3),trace(:,3),'marker','*');
    hold(SingleTraceH,'off');
    xlim(SingleTraceH,[0,(max(trace(:,1))-1)*I.xyt_dim(3)]);
    box(SingleTraceH,'off');
    % axis(SingleTraceH,'tight'); 
    grid(SingleTraceH,'on');
    title(SingleTraceH,'Spark Trace');
    xlabel(SingleTraceH,'Time (ms)'); ylabel(SingleTraceH,'Fluorescence (a.u.)');
    
    %% Full trace.
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        FullTraceH=subplot(2,20,21:34,'Parent',Live_fig);
    end
    % Plot full pixel transient.
    plot(FullTraceH,(0:(numel(tracefull)-1))*I.xyt_dim(3)/1000,tracefull,...
        'DisplayName','Denoised','Color',[0,0,1]);
    hold(FullTraceH,'all')
    
    % Plot spark location/postion.
    plot(FullTraceH,ones(1,2)*t*I.xyt_dim(3)/1000,[sparkBgr*0.9, (sparkBgr+sparkAmp)*1.1],...
        'DisplayName','Location','Color',[1,0,1]);
    
    % Plot a box to hold the current spark.
    MarkerY=(round([(t1-I.SparkLim(1)/I.xyt_dim(3))*I.xyt_dim(3),...
        (t2+I.SparkLim(2)/I.xyt_dim(3))*I.xyt_dim(3)])); MarkerY=MarkerY/1000;
    plot(FullTraceH,[MarkerY(1),MarkerY(2),MarkerY(2),MarkerY(1),MarkerY(1)],...
        [sparkBgr,sparkBgr,sparkAmp+sparkBgr,sparkAmp+sparkBgr,sparkBgr],...
        'Color',[1,0,0]);
    
    xlim_max=ceil((numel(tracefull)-1)*I.xyt_dim(3)/100)/10;
    xlim(FullTraceH,[0 xlim_max]);
    set(FullTraceH,'XTick',0:ceil(xlim_max/10):xlim_max)
    plotmax=max([max(tracefull),sparkAmp+sparkBgr,sparkBgr]);
    plotmin=min([min(tracefull),sparkAmp+sparkBgr,sparkBgr]);
    plotmin=round(plotmin-abs(plotmin)*0.02); plotmax=round(plotmax+abs(plotmax)*0.02);
    ylim(FullTraceH,[plotmin plotmax]);
    set(FullTraceH,'YTick',round(plotmin:ceil((plotmax-plotmin)/5):plotmax))
    hold(FullTraceH,'off');
    xlabel(FullTraceH,'Time (s)');ylabel(FullTraceH,'Fluorescence (a.u.)');
    box(FullTraceH,'off');
    grid(FullTraceH,'on');
    title(FullTraceH,['Spark No.' num2str(ID) '    |    Central pixel transient: dAmp=' ...
        num2str(round(sparkAmp)) '    |    Bgr=' num2str(round(sparkBgr)) ...
        '    |    dF/F0=' num2str(round(sparkFF0*100)/100)]);
    if ishandle(Live_fig)==0
        StatusGood='Stop';
        return;
    else
        drawnow update;
    end

    %% Manual Check.
    if ManualCheck
        if isempty(MFquestdlgPos)
            MFquestdlgPos=dlgSize;
        end
        [status,MFquestdlgPos]=MFquestdlg(MFquestdlgPos,['Checking spark ' num2str(ID) ...
            ': Accept it as a good spark?'],'Single Spark Checking',...
            'Accept','Reject','Stop','Accept');
        StatusGood=status;
        if strcmpi(status,'Stop')
            return;
        end

        msg=['    No.%g  I(%0.2f,%0.2f,%0.2f)  '...         %No., Center
        '%4.2f  %0.4f  %0.4f  %0.4f  ' ...                  %Spark Mass, FWHM, FWHM_max and FWHM_min.
        '%1.4f  %1.4f  %1.4f  %4.2f  %4.2f  %1.2f %s\n'];   %AmpFF0,Tau,Sigma,AmpInt,BgrInt,SparkToMembraneDistance
        fprintf(msg,ID,I.SparkPos(ID,7),I.SparkPos(ID,8),I.SparkPos(ID,9),...
            I.SparkProperty(ID,10),I.SparkProperty(ID,1),I.SparkProperty(ID,2),I.SparkProperty(ID,3),...
            I.SparkProperty(ID,6),I.SparkProperty(ID,7),I.SparkProperty(ID,8),...
            I.SparkProperty(ID,4),I.SparkProperty(ID,5),I.SparkProperty(ID,11),status);
    end
    
end










%%
function Pos=SparkFigSize
    Screen=get(0,'ScreenSize');
    FigHeight=round(Screen(4)*0.7);
    FigWidth=round(2.2857*FigHeight);
    if FigWidth>Screen(3)*0.95
        FigWidth=round(Screen(3)*0.95);
        FigHeight=round(FigWidth/2.2857);
    end
    Pos=[round((Screen(3)-FigWidth)*0.5) round((Screen(4)-FigHeight)*0.5) FigWidth FigHeight];
end

function Pos=dlgSize
    Screen=get(0,'ScreenSize');
    FigHeight=round(Screen(4)*0.7);
    
    Pos=[(Screen(3)-267)/Screen(3)/4 ((round((Screen(4)-FigHeight)*0.7)-67)*0.3)/Screen(4)];
end


function BlackToRed=coolcolor
    BlackToRed=[0 0 0;0 0 0.0666666701436043;0 0 0.133333340287209;0 0 0.200000002980232;...
        0 0 0.266666680574417;0 0 0.333333343267441;0 0 0.400000005960464;0 0 0.466666668653488;...
        0 0 0.533333361148834;0 0 0.600000023841858;0 0 0.666666686534882;0 0 0.733333349227905;...
        0 0 0.800000011920929;0 0 0.866666674613953;0 0 0.933333337306976;0 0 1;0 0.0625 1;...
        0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;...
        0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;...
        0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;...
        0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;...
        0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;...
        1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;...
        1 0.0625 0;1 0 0];
end

function [bigImage,bw] = montage_Modified(I,mSize)
    %   Copyright 1993-2011 The MathWorks, Inc.
    %   $Revision: 1.1.8.17 $  $Date: 2011/10/11 15:49:16 $
    
    indices=1:size(I,4);
    nFrames = numel(indices);
    nRows = size(I,1);
    nCols = size(I,2);
    
    montageSize = calculateMontageSize(mSize);
    
    bigImage = createMontageImage;
    bigImage = bigImage(1:(size(bigImage,1)-1),1:(size(bigImage,2)-1));
    bw=false(size(bigImage,1),size(bigImage,2));
    bw((nRows+1):(1+nRows):end,:)=true;
    bw(:,(nCols+1):(1+nCols):end)=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function montageSize = calculateMontageSize(mSize)
        montageSize = mSize;
        nanIdx = isnan(mSize);
        montageSize(nanIdx) = ceil(nFrames / mSize(~nanIdx));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function b = createMontageImage
        nMontageRows = montageSize(1);
        nMontageCols = montageSize(2);
        nBands = size(I, 3);
        
        sizeOfBigImage = [nMontageRows*(nRows+1) nMontageCols*(nCols+1) nBands 1];
        if islogical(I)
            b = false(sizeOfBigImage);
        else
            b = zeros(sizeOfBigImage, class(I));
        end
        
        rows = 1 : nRows;
        cols = 1 : nCols;
        k = 1;
        
        for j = 0 : nMontageCols-1
            for i = 0 : nMontageRows-1
                if k <= nFrames
                    b(rows + i * (1+nRows), cols + j * (1+nCols), :) = I(:,:,:,indices(k));
                else
                    return;
                end
                k = k + 1;
            end
        end
    end
end