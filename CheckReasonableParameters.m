function [SparkPosKept,SparkPropertyKept,SparkPosRejected,SparkPropertyRejected]=...
        CheckReasonableParameters(SparkPos,SparkProperty,GauAmpLim,TauLim,FWHMLim)
    %% Check Reasonable Parameters.

    % ======================================
    % Just for internal usage. Qinghai Tian.
    % ======================================
    %
    % Input:
    %  SparkPos=[x1,x2,y1,y2,z1,z2,xc,yc,z_onset,No.];
    %  SparkProperty_oldversion=[FWHM,FWHM_R2,dAmp,Bgr,Tau,Sigma,R2,dAmp_dF/F0,DetctionMass,dF/F0_Mass, No.];
    %  SparkProperty=[FWHM,MajorAxisLength,MinorAxisLength,Amp,Bgr,Amp/Bgr,Tau,Sigma,DetctionMass,dF/F0_Mass, No.];
    %                   1        2              3           4   5     6     7    8        9         10        11 
 
    %% Checking fitting results.
    fprintf('  0%%');
    TotalNum=size(SparkProperty,1);
    
    SparkPropertyKept=zeros(TotalNum,size(SparkProperty,2));
    SparkPropertyRejected=zeros(TotalNum,size(SparkProperty,2)+3);
    SparkPosKept=zeros(size(SparkPos));
    SparkPosRejected=zeros(size(SparkPos));

    
    KeptNo=0;RejectedNo=0;
    for k=1:TotalNum
        if SparkProperty(k,3)~=0
            AmpCheck=((SparkProperty(k,6) <= GauAmpLim) & (SparkProperty(k,6) >0));   % Amplitude check
            FWHMCheck=(SparkProperty(k,1)>FWHMLim(1)) & (SparkProperty(k,1)<FWHMLim(2)); % FWHM check
            TauCheck=(SparkProperty(k,7)>0 & SparkProperty(k,7)<TauLim);
            
            if AmpCheck && FWHMCheck && TauCheck 
                KeptNo=KeptNo+1;
                SparkPropertyKept(KeptNo,:)=SparkProperty(k,:);
                SparkPosKept(KeptNo,:)=SparkPos(k,:);

            else
                RejectedNo=RejectedNo+1;
                SparkPropertyRejected(RejectedNo,:)=...
                    [SparkProperty(k,:) AmpCheck FWHMCheck TauCheck];
                SparkPosRejected(RejectedNo,:)=SparkPos(k,:);
            end
        else
        end
        fprintf('\b\b\b\b%3.0f%%',floor(k/TotalNum*100));
    end
    SparkPosKept=SparkPosKept(1:KeptNo,:);
    SparkPropertyKept=SparkPropertyKept(1:KeptNo,:);
    
    SparkPosRejected=SparkPosRejected(1:RejectedNo,:);
    SparkPropertyRejected=SparkPropertyRejected(1:RejectedNo,:);
    fprintf('. %g final sparks',KeptNo);
end
