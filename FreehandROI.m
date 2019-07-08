function  FinalBW = FreehandROI(Iall,mask,strTitle)
    % FreehandROI allows you to draw multiple ROIs from an input image.
    if size(Iall,3)>1; Iall=mean(Iall,3); end
    if ~isfloat(Iall); Iall=single(Iall); end
    
    Iallsort=sort(Iall(:));
    FigHandle=figure;
    try
        currFigContentHandle=imagesc(Iall,[Iallsort(round(numel(Iallsort)*0.01)) Iallsort(round(numel(Iallsort)*0.99))]);
    catch
        currFigContentHandle=imagesc(Iall);
    end
    currentFigAx=get(gcf,'CurrentAxes');
    axis(currentFigAx,'image','off');
    if nargin==3
        title(currentFigAx,strTitle,'Interpreter','none');
    else
        title(currentFigAx,'ROI Drawing');
    end
    set(FigHandle,'Position',FigSize(size(Iall,1),size(Iall,2)));
    colormap(currentFigAx,coolcolor);
    
    [m,n]=size(Iall);
    
    if (nargin>=2)&&(size(Iall,1)==size(mask,1))&&(size(Iall,2)==size(mask,2))
        FinalBW=mask;
    else
        figure(FigHandle);
        h = imfreehand;
        pos = getPosition(h);
        FinalBW=poly2mask(pos(:,1), pos(:,2), m, n);
        delete(h);
    end
    
    CellEdge=FinalBW-ordfilt2(FinalBW, 1, [0 1 0;1 1 1;0 1 0],'zeros')==1;
    IallBW=Iall;IallBW(CellEdge)=Iallsort(round(numel(Iallsort)*0.99));
    set(currFigContentHandle,'CData',IallBW);
    
    
    bwContinue=true;
    msgPos=dlgSize;
    while bwContinue
        figure(FigHandle);
        [status,msgPos]=MFquestdlg(msgPos,sprintf(['Are you satisfied with this mask?' ...
            ' "Yes" to finish drawing.\n \nButton "Add" to add a region, ' ...
            '"Remove" to remove a region.']),'Decision','Yes','Add','Remove','Yes');
        if strcmpi(status,'Yes')
            break
        else
            h = imfreehand;
            pos = getPosition(h);
            BW = poly2mask(pos(:,1), pos(:,2), m, n);
        end
        
        if strcmpi(status,'Add')
            FinalBW = FinalBW | BW;
        elseif strcmpi(status,'Remove')
            FinalBW = FinalBW & (~BW);
        end
        CellEdge=FinalBW-ordfilt2(FinalBW, 1, [0 1 0;1 1 1;0 1 0],'zeros')==1;
        IallBW=Iall;IallBW(CellEdge)=Iallsort(round(numel(Iallsort)*0.99));
        set(currFigContentHandle,'CData',IallBW);
        delete(h);drawnow update
    end
    close(FigHandle);drawnow;
end


function Pos=FigSize(x,y)
    
    Screen=get(0,'ScreenSize');
    
    if Screen(3)/y < Screen(4)/x
        FigWidth=Screen(3)*0.9;
        FigHeight=FigWidth/y*x;
    else
        FigHeight=Screen(4)*0.8;
        FigWidth=FigHeight/x*y;
    end
    
    Pos=[round((Screen(3)-FigWidth)*0.5) round((Screen(4)-FigHeight)*0.5) FigWidth FigHeight];
end


function Pos=dlgSize
    Screen=get(0,'ScreenSize');
    FigHeight=round(Screen(4)*0.7);
    
    Pos=[(Screen(3)-267)/Screen(3)/4 ((round((Screen(4)-FigHeight)*0.7)-67)*0.3)/Screen(4)];
end


function BlackToRed=coolcolor
    BlackToRed=[0 0 1;0 0.0476190485060215 1;0 0.095238097012043 1;...
        0 0.142857149243355 1;0 0.190476194024086 1;0 0.238095238804817 1;...
        0 0.28571429848671 1;0 0.333333343267441 1;0 0.380952388048172 1;...
        0 0.428571432828903 1;0 0.476190477609634 1;0 0.523809552192688 1;...
        0 0.571428596973419 1;0 0.61904764175415 1;0 0.666666686534882 1;...
        0 0.714285731315613 1;0 0.761904776096344 1;0 0.809523820877075 1;...
        0 0.857142865657806 1;0 0.904761910438538 1;0 0.952380955219269 1;...
        0 1 1;0.0476190485060215 1 0.952380955219269;...
        0.095238097012043 1 0.904761910438538;...
        0.142857149243355 1 0.857142865657806;...
        0.190476194024086 1 0.809523820877075;...
        0.238095238804817 1 0.761904776096344;...
        0.28571429848671 1 0.714285731315613;...
        0.333333343267441 1 0.666666686534882;...
        0.380952388048172 1 0.61904764175415;...
        0.428571432828903 1 0.571428596973419;...
        0.476190477609634 1 0.523809552192688;...
        0.523809552192688 1 0.476190477609634;...
        0.571428596973419 1 0.428571432828903;...
        0.61904764175415 1 0.380952388048172;...
        0.666666686534882 1 0.333333343267441;...
        0.714285731315613 1 0.28571429848671;...
        0.761904776096344 1 0.238095238804817;...
        0.809523820877075 1 0.190476194024086;...
        0.857142865657806 1 0.142857149243355;...
        0.904761910438538 1 0.095238097012043;...
        0.952380955219269 1 0.0476190485060215;...
        1 1 0;1 0.952380955219269 0;1 0.904761910438538 0;...
        1 0.857142865657806 0;1 0.809523820877075 0;...
        1 0.761904776096344 0;1 0.714285731315613 0;...
        1 0.666666686534882 0;1 0.61904764175415 0;...
        1 0.571428596973419 0;1 0.523809552192688 0;...
        1 0.476190477609634 0;1 0.428571432828903 0;...
        1 0.380952388048172 0;1 0.333333343267441 0;...
        1 0.28571429848671 0;1 0.238095238804817 0;...
        1 0.190476194024086 0;1 0.142857149243355 0;...
        1 0.095238097012043 0;1 0.0476190485060215 0;1 0 0];
end