function Threshold=AdaptiveThreshold(DataCV,denoise)
    switch denoise
        case 1
            % Without imdilate([0,1,0;1,1,1;0,1,0]),imerode and imfill('holes')
            % Puredenoise on xy, xt and yt
            x=[0.0280350595712662;0.0274587161839008;0.0252197682857513;...
                0.0216009318828583;0.0181888807564974;0.0145628452301025;...
                0.0127790514379740;0.0104524996131659;0.00820331741124392;...
                0.00674492213875055;0.00576218916103244;0.00491002667695284;...
                0.00428543100133538;0.00351920211687684;0.00268866331316531;...
                0.00188843207433820;0.00132414000108838;0.00106183555908501;...
                0.000807940494269133;0.000598307175096124];
            y=[7.93663332817856;7.60233909124705;7.62693846750941;8.06295151655582;...
                8.82085761782565;8.94049778853035;9.14535653396641;8.92516435472977;...
                8.73925658013017;8.16834408459598;8.43713108948362;8.12322224082031;...
                8.37513621775521;8.02790814143723;8.14564069488819;8.46669962057445;...
                8.41119892209994;8.66307573902548;8.73332472949136;9.07606684547858];
            
            % PUREdenoise xy only.
            % x=[0.133524239063263;0.113800361752510;0.100217044353485;0.0840259045362473;...
            %     0.0659703314304352;0.0525200106203556;0.0452923029661179;0.0390192121267319;...
            %     0.0341418758034706;0.0282632783055305;0.0222283508628607;0.0162998139858246;...
            %     0.0123412311077118;0.0106835132464767;0.00905625335872173;0.00755277508869767];
            % y=[5.13510091987235;4.95932055150989;4.98964345850957;4.97283450511497;...
            %     4.92505295135496;4.86133470881836;4.84417317719228;4.81891158178368;...
            %     4.77650874249642;4.78276556161098;4.70460788555079;4.69597225855308;...
            %     4.64584347756503;4.59166848135837;4.51723176085423;4.40669040197122];
            
        otherwise % CANDLEdenoise or no denoise.
            x=[0.0325001068413258;0.0242019146680832;0.0188944879919291;0.0141100250184536;...
                0.0109231835231185;0.00921791139990091;0.00786810275167227;0.00678965914994478;...
                0.00553421862423420;0.00429173326119781;0.00304904114454985;0.00217796978540719];
            y=[6.67768360588211;7.15140882419229;6.76013440148786;6.58218105498078;6.53767193703302;...
                6.57390704662395;6.52718443255285;6.53752273110893;6.41318941916338;6.50440019669000;...
                6.47764104644472;6.20828271628213];
    end
    x_min=min(x);
    x_max=max(x);
    DataCV(DataCV<x_min)=x_min;
    DataCV(DataCV>x_max)=x_max;
    
    Threshold = interp1(x,y,DataCV,'linear','extrap');
end

% function Threshold=AdaptiveThreshold(DataCV) % No imdilate, imerode and imfill
%     x=[0.0280884355306625;0.0267239101231098;0.0251312311738729;0.0221580881625414;...
%         0.0178016899153590;0.0147603303194046;0.0124648613855243;0.0103105427697301;...
%         0.00817436818033457;0.00664802500978112;0.00569266080856323;0.00496280053630471;...
%         0.00428770272992551;0.00346726621501148;0.00269968179054558;0.00190746749285609;...
%         0.00134606217034161;0.00107926520286128;0.000821087247459218;0.000591837975662202];
%     y=[7.92258961314196;7.48153602954865;7.88113231322362;7.79597594183423;...
%         8.87478162618359;8.61371705897528;8.92629180364144;9.17064534077250;...
%         9.01634951167465;8.69751739174959;8.48150733733186;8.28129477963109;...
%         7.98671040449138;8.08855525202311;8.03066084064988;8.16014396967764;...
%         8.17831797234835;8.45538329351248;8.40431267890245;9.11157539053375];
%
%     x_min=min(x);
%     x_max=max(x);
%     DataCV(DataCV<x_min)=x_min;
%     DataCV(DataCV>x_max)=x_max;
%
%     Threshold = interp1(x,y,DataCV,'linear','extrap');
% end



% function Threshold=AdaptiveThreshold(PhotonNum) % No imdilate, imerode and imfill
%     x=[0.5;0.6;0.7;0.8;0.9;1;1.5;2;3;5;8;11.;15;20;30;50;100;200;300;500;1000];
%     y=[12.6337965478292;8.70770211162130;8.40677554376737;8.30105426259988;...
%     9.08552021140613;10.2571329599866;10.7656268946016;10.8054367376458;11.1135846697754;...
%     10.9567053115862;10.4157247126620;10.3558305259946;9.96422223710604;9.79377607212849;...
%     9.97714634412640;9.91222739044396;10.2425298604482;10.2849126963120;10.8345119354181;...
%     11.0484783795897;12.0530367167779];
%     x_min=min(x);
%     x_max=max(x);
%     PhotonNum(PhotonNum<x_min)=x_min;
%     PhotonNum(PhotonNum>x_max)=x_max;
%
%     Threshold = interp1(x,y,PhotonNum,'linear','extrap');
% end

% PUREdenoise on xy, xt, and yt, respectively.
% function Threshold=AdaptiveThreshold(DataCV)
%     TH=@(Y0,Plateau,K,X) Y0 + (1-exp(-K*X)) * (Plateau - Y0);
%
%     Y0=-7.203;
%     Plateau=11.44;
%     K=1940;
%
%     Threshold=TH(Y0,Plateau,K,DataCV);
%
%     Threshold(isnan(Threshold))=Plateau;
%     Threshold(DataCV<0)=Plateau;
% end


% function Threshold=AdaptiveThreshold(DataCV)  % Sinlge frame based
% puredenoise.
%     TH=@(Y0,Plateau,K,X) Y0 + (1-exp(-K*X)) * (Plateau - Y0);
%
%     Y0=-3.236;
%     Plateau=5.137;
%     K=103.3;
%
%     Threshold=TH(Y0,Plateau,K,DataCV);
%
%     Threshold(isnan(Threshold))=Plateau;
%     Threshold(DataCV<0)=Plateau;
% end