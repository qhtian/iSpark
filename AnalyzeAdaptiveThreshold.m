function finalTH=AnalyzeAdaptiveThreshold(TH)
    % AnalyzeAdaptiveThreshold calcuates the suitble threshold from simulation.
    %
    % Syntax:
    %   finalTH=AnalyzeAdaptiveThreshold(SavedFinalThNum);
    %
    % finalTH should have the following format:
    %   FinalTH(1,:): input gain value;
    %   FinalTH(2,:): raw DataCV;
    %   FinalTH(3,:): detected gain value;
    %   FinalTH(4,:): denoised DataCV;
    %   FinalTH(5,:): optimal threshold;
    %   -----------------------------
    % SavedFinalThNum can the output directly from SimCalAdpTH.
    %
    % See also SimCalAdpTH.
    
    
    numOfCalculation=size(TH.SimCalResults,2)/2;
    finalTH=zeros(numOfCalculation,5);
    
    for k=1:numOfCalculation
        currFalsePositive=TH.SimCalResults(5:end,k*2);
        currThreshold=TH.SimCalResults(5:end,k*2-1);
        bw=currFalsePositive>5;
        sw=false;
        for j=2:numel(bw)
            if sw
                bw(j)=false;
                continue;
            end
            if currFalsePositive(j)<currFalsePositive(j-1)
                bw(j)=false;
                sw=true;
            end
        end
        clear('sw')
        currFalsePositive=log2(currFalsePositive(bw));
        currThreshold=currThreshold(bw);
        
        clear('bw')
        try
            p=robustfit(currThreshold,currFalsePositive);
        catch
            p=[0.001,0.001];
        end
        h=plot(currThreshold,currFalsePositive,'o');drawnow;pause(0.1)
        hold('all');
        plot(currThreshold,p(1)+p(2)*currThreshold,'color',get(h,'color'));
        % xlim([0,8])
        finalTH(k,1)=TH.SimCalResults(1,k*2);
        finalTH(k,2)=TH.SimCalResults(2,k*2);
        finalTH(k,3)=TH.SimCalResults(3,k*2);
        finalTH(k,4)=TH.SimCalResults(4,k*2);
        finalTH(k,5)=-p(1)/p(2);
    end
    xlabel('Threshold');
    ylabel('False Positive Events');
end
