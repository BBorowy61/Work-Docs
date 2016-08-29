% McHenry PJM Analysis
% 
% Bogdan Borowy
% EDF-RE, 2016

tstart = tic;

Get_Data;

%%%%%%%%  User Entry Data %%%%%%%%%%%%
RR=EnergyRRatAGCMWMWmin/6;

UnitAREG = 20*ones(size(RR));
FleetTREG = 2000*ones(size(RR));

UseLoadBP = false;

SlopeToUse = 'Percent';
% SlopeToUse = 'Signal'

HH2 = 0;  % unknown condition in $HH$2
HH5 = 0;  % unknown condition in $HH$5

SlopeLimitVal = 0;
StdDevLimitVal = 0.05;

A=0.3333;
B=0.3333;
C=0.3333;
%%%%%%%%  &&&&&&&&&&&&&&&&&&&&& %%%%%%%%%%%%

%[G]
if(UseLoadBP == true)|(abs(RawECOBP-Vshift(RampedEcoBP,-1))<RR),
    RampedEcoBP = RawECOBP;
else
    RampedEcoBP = Vshift(RampedEcoBP,-1)+sign(RawECOBP-Vshift(RampedEcoBP,-1)).*RR;
end

%[K]
RawUnitReg = zeros(size(RR));
if(FleetTREG~=0),
    RawUnitReg = UnitAREG .* FleetRegulationControlSignal./FleetTREG;
end

%[L]
RawEcoBP_RawReg = RawECOBP + RawUnitReg;

%[M]
UnitResp = AGCMW - RampedEcoBP;

%[N]

%[AD - BH]
%Correlation Signals to Response.
% Every 5 minutes we calculate new correlation coefficient.  Ts = 10s. Then
% number of samples to take for correlation calcs is N = 5 min * 60 sec /Ts
N = 30;
Ts=10;
H1int = 60 * 60 / Ts;  %analysis interval is 1 hour
%Correlation matrix.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the commanded signal and the response - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
CorrelSig2Respsec = zeros(H1int,M);  

for j=1:M,
    S = j - 1;  % 10 sec time shifts up to 300 sec (5 min) shift
     
    temp = Vshift(UnitResp,S);
    for i = 1:H1int,
        temp1 = Vshift(RawUnitReg,i-1);
        temp2 = Vshift(temp,i-1);

        CorrelSig2RespsecX = corrcoef(temp1(1:N+1), temp2(1:N+1));
        CorrelSig2Respsec(i,j) = CorrelSig2RespsecX(1,2);
    end
end
%Plot the Correlation Matrix
[X,Y] = meshgrid(1:10:3600,0:10:300);

figure
subplot(2,1,1);
surf(X,Y,CorrelSig2Respsec');
axis([0 3600 0 300 -1 1]);
title('Correlatino Matrix');
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap copper;
colorbar('westoutside');
%colorbar('Ticks',[-0.9,-0.5,0,0.5,0.9],...
%         'TickLabels',{'Bad','Poor','Acceptable','OK','Good'})

subplot(2,1,2);
contour(X,Y,CorrelSig2Respsec');
axis([0 3600 0 300]);
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap copper;
colorbar('westoutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [BJ]
temp = zeros(size(FleetRegulationControlSignal));
if(FleetTREG ~= 0),
    temp = FleetRegulationControlSignal./FleetTREG;
end

RegSignalPercent = min(max(temp,-1),1);  %make sure the signal is within <-1;1> domain
% [BK]
% slope(RawUnitReg,Interval) = cov(Interval,RawUnitReg)/var(Interval);
Interval = (1:length(RR))';

FiveMinuteSignalSlope = zeros(H1int,1);
for i = 1:H1int,
    temp1 = Vshift(RawUnitReg,i-1);
    temp2 = Vshift(Interval,i-1);
    slope= cov(temp1(1:N+1),temp2(1:N+1))/var(temp2(1:N+1));
    FiveMinuteSignalSlope(i) = slope(1,2);
end

% [BL]
% slope(RegSignalPercent,Interval) = cov(Interval,RegSignalPercent)/var(Interval);
Interval = (1:length(RR))';
FiveMinuteSignalPercentSlope = zeros(H1int,1);
for i = 1:H1int,
    temp1 = Vshift(RegSignalPercent,i-1);
    temp2 = Vshift(Interval,i-1);
    slope= cov(temp1(1:N+1),temp2(1:N+1))/var(temp2(1:N+1));
    FiveMinuteSignalPercentSlope(i) = slope(1,2);
end

%[BM] SlopeUsed
if(strcmpi(SlopeToUse,'SIGNAL')),
    SlopeUsed = FiveMinuteSignalSlope;
else
    SlopeUsed = FiveMinuteSignalPercentSlope;
end
    
% [BN] SlopeLimit
SlopeLimit = SlopeLimitVal*ones(size(FiveMinuteSignalSlope));
% [BO] -SlopeLimit
SlopeLimitNeg = -SlopeLimit;

% [BP] StdDevUsed
StdDevUsed = zeros(H1int,1);
for i = 1:H1int,
    temp1 = Vshift(RegSignalPercent,i-1);
    %slope= std(temp1(1:N+1));
    StdDevUsed(i) = std(temp1(1:N+1));%slope(1,2);
end

% [BQ] StdDevLimit
StdDevLimit = StdDevLimitVal*ones(size(FiveMinuteSignalSlope));

% [BS]-[CW]  Response slopes
%   Slope matrix.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the UnitResp and the Measurement Interval - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
ResponseSlope = zeros(H1int,M);  
Interval = (1:length(RR))';
for j=1:M,
    S = j - 1;  % 10 sec time shifts up to 300 sec (5 min) shift 
    temp = Vshift(UnitResp,S);
    for i = 1:H1int,
        tempy = Vshift(temp,i-1);
        tempx = Vshift(Interval,i-1);
        slope= cov(Interval(1:N+1),tempy(1:N+1))/var(Interval(1:N+1));
        ResponseSlope(i,j) = slope(1,2);
    end
end
%Plot the Correlation Matrix
[X,Y] = meshgrid(1:10:3600,0:10:300);

figure
subplot(2,1,1);
surf(X,Y,ResponseSlope');
axis([0 3600 0 300 -1 1]);
title('ResponseSlope Matrix');
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap jet;
colorbar('westoutside');
%colorbar('Ticks',[-0.9,-0.5,0,0.5,0.9],...
%         'TickLabels',{'Bad','Poor','Acceptable','OK','Good'})

subplot(2,1,2);
contour(X,Y,ResponseSlope');
axis([0 3600 0 300]);
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap copper;
colorbar('westoutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CY] - [EC] Accuracy Scoring Mode Matrix
%   Scoring matrix.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the UnitResp and the Measurement Interval - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
AccuracyScoringMode = cell(H1int,M);  
for j=1:M,
    for i = 1:H1int,
        if((abs(ResponseSlope(i,j)) < SlopeLimitVal)||...
                (StdDevUsed(i) < StdDevLimitVal)),
            AccuracyScoringMode{i,j} = 'Slope';
        else
            AccuracyScoringMode{i,j} = 'Correl';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [EE] - [FE] Accuracy Score Matrix
%   Scoring matrix.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the UnitResp and the Measurement Interval - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
AccuracyScore = zeros(H1int,M);  
for j=1:M,
    for i = 1:H1int,
        if(strcmp(AccuracyScoringMode{i,j},'Correl')),
            temp = CorrelSig2Respsec(i,j);
        else
            temp = 1 - abs(SlopeUsed(i)-ResponseSlope(i,j));
        end
        AccuracyScore(i,j) = max(min(temp,1),0);
    end
end
%Plot the Correlation Matrix
[X,Y] = meshgrid(1:10:3600,0:10:300);

figure
subplot(2,1,1);
surf(X,Y,AccuracyScore');
axis([0 3600 0 300 -1 1]);
title('AccuracyScore Matrix');
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap jet;
colorbar('westoutside');
%colorbar('Ticks',[-0.9,-0.5,0,0.5,0.9],...
%         'TickLabels',{'Bad','Poor','Acceptable','OK','Good'})

subplot(2,1,2);
contour(X,Y,AccuracyScore');
axis([0 3600 0 300]);
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap copper;
colorbar('westoutside');

% [FK] - [FM] Max Accuracy Score, Offset, Accuracy Method Used
MaxAccuracyScore = zeros(H1int,1);
Offset = zeros(H1int,1);
AccuracyMethodUsed = cell(H1int,1);
for i = 1:H1int,
    if FleetTREG(i)~=0,
        MaxAccuracyScore(i) = max(AccuracyScore(i,:));
        Offset(i) = find(MaxAccuracyScore(i)==AccuracyScore(i,:))-1;
        AccuracyMethodUsed{i} = AccuracyScoringMode{Offset(i)+1};
    end
end

% [FN]  Delay score
DelayScore = zeros(H1int,1);
for i = 1:H1int,
    if MaxAccuracyScore(i)~=0,
        DelayScore(i) = min((1-(Offset(i)-1)/30),1);
    end
end


% [FP] - [GT]
%   Weighted Acc & Delay.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the UnitResp and the Measurement Interval - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
WeightedAccANDDelay = zeros(H1int,M);  
for j=1:M,
    for i = 1:H1int,
       WeightedAccANDDelay(i,j) = A * AccuracyScore(i,j) + B *  min((1-(j-2)/30),1);
    end
end
%Plot the Correlation Matrix
[X,Y] = meshgrid(1:10:3600,0:10:300);

figure
subplot(2,1,1);
surf(X,Y,WeightedAccANDDelay');
axis([0 3600 0 300 -1 1]);
title('WeightedAccANDDelay Matrix');
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap jet;
colorbar('westoutside');
%colorbar('Ticks',[-0.9,-0.5,0,0.5,0.9],...
%         'TickLabels',{'Bad','Poor','Acceptable','OK','Good'})

subplot(2,1,2);
contour(X,Y,WeightedAccANDDelay');
axis([0 3600 0 300]);
xlabel('10 sec moving window');
ylabel('10 set time shift');
colormap copper;
colorbar('westoutside');

% [GV] - [GW]  MaxWeightedAccDel + Offset1 
MaxWeightedAccDel = zeros(H1int,1);
Offset1 = zeros(H1int,1);
for i = 1:H1int,
    MaxWeightedAccDel(i) = max(WeightedAccANDDelay(i,:));
    Offset1(i) = find(MaxWeightedAccDel(i)==WeightedAccANDDelay(i,:))-1;
end


% [GX]  AccuracyScoreAdjust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [GY] AccuracyScoreAdjust Matrix
%   Scoring matrix.
%   1) Each Column number represents a multiple of 10 sec time shift between
%   the UnitResp and the Measurement Interval - from 0 up to 5 minutes relative
%   shifts.
% 
%   2) Each Row represents 5 minutes moving window in 10 sec intervals of both
%   response and the commanded singnals.
M = 5*60/10+1;  % number of columns of the Correlatino Matrix
AccuracyScoreAdjust = zeros(H1int,1); 

for i = 1:H1int,
    AccuracyScoreAdjust(i) = AccuracyScore(i,Offset1(i)+1);
    temp1 = Vshift(UnitResp,i-1);
    temp2 = Vshift(UnitAREG,i-1);
    if((max(temp1(1:2*N+1)) == min(temp1(1:2*N+1))) || ...
            (AccuracyScore(i,Offset1(i)+1)<0.0000001)||...
            (FleetTREG(i)==0)),
        AccuracyScoreAdjust(i)=0;
    elseif (min(temp2(1:N+1))==0)
          AccuracyScoreAdjust(i)=NaN;
    end    
end

% [GY]  AccuracyMethodUsedAdj

AccuracyMethodUsedAdj = cell(H1int,1);
for i = 1:H1int,
    AccuracyMethodUsedAdj{i} = AccuracyScoringMode{Offset1(i)+1};
end


% [GZ]  DelayScoreAdj
DelayScoreAdj = zeros(H1int,1); 
for i = 1:H1int,
    if min(UnitAREG)==0,
        DelayScoreAdj(i) = NaN;
    elseif AccuracyScoreAdjust(i)~=0,
        DelayScoreAdj(i) = min((1-(Offset1(i)-1)/30),1);
    end
end

%[HB]  AbsoluteSignal
AbsoluteSignal = abs(RawUnitReg);

% [HC] AverageAbsSignal
AverageAbsSignal = mean(AbsoluteSignal);

% [HE] PrecisScoreAREG
PrecisScoreAREG = zeros(H1int,1); 
temp = Vshift(UnitResp,1);
for i = 1:H1int,
    if FleetTREG(i) ~= 0,
        PrecisScoreAREG(i) = min(max(1-abs(temp(i) - RawUnitReg(i))/UnitAREG(i),0),1);
    end
end

%[HF]  PrecisScoreAbsSig

PrecisScoreAbsSig = zeros(H1int,1); 
temp = Vshift(UnitResp,1);
if (max(UnitResp) ~= min(UnitResp)),
    if(min(UnitAREG) == 0),
        PrecisScoreAbsSig = ones(size(PrecisScoreAREG)) * NaN;
    else
        for i = 1:H1int,
            if FleetTREG(i) ~= 0,
                PrecisScoreAbsSig(i) = min(max(1-abs(temp(i)-RawUnitReg(i))/AverageAbsSignal,0),1); 
            end
        end
    end
end

% [HH]  AccuracyAdjust2
        
if HH2 == 1,
    AccuracyAdjust2 = MaxAccuracyScore;
else
    AccuracyAdjust2 = AccuracyScoreAdjust;
end

% [HI] AccuracyMethodUsedAdj2
if HH2 == 1,
    AccuracyMethodUsedAdj2 = AccuracyMethodUsed;
else
    AccuracyMethodUsedAdj2 = AccuracyMethodUsedAdj;
end

% [HJ]  AccuracyMethodNumber
AccuracyMethodNumber = zeros(H1int,1);
for i = 1:H1int,
    if(strcmp(AccuracyMethodUsedAdj2{i},'Correl')),
        AccuracyMethodNumber(i) = 1;
    elseif(strcmp(AccuracyMethodUsedAdj2{i},'Slope')),
        AccuracyMethodNumber(i) = 2;
    end
end

% [HK]  DelayScoreAdj2
if HH2 == 1,
    DelayScoreAdj2 = DelayScore;
else
    DelayScoreAdj2 = DelayScoreAdj;
end

% [HL]  PrecisionScore
PrecisionScore = PrecisScoreAbsSig;
if HH5 == 1,
    PrecisionScore = PrecisScoreAREG;
end

% [HP] Mileage
Mileage = zeros(H1int,1);
temp = Vshift(RegSignalPercent,-1);
for i = 1:H1int,
    Mileage(i) = abs(RegSignalPercent(i) - temp(i));
end

% [HR] StdDevofSlopes
StdDevofSlopes = zeros(H1int,1);

for i = 1:H1int,
    temp1 = Vshift(SlopeUsed,i-1);
    stdsslp= std(temp1(1:N+1));
    StdDevofSlopes(i) = stdsslp;
end

% [HS] StdDevOfPercentSignal
StdDevOfPercentSignal = zeros(H1int,1);

for i = 1:H1int,
    temp1 = Vshift(RegSignalPercent,i-1);
    stdsslp= std(temp1(1:N+1));
    StdDevOfPercentSignal(i) = stdsslp;
end

telapsed = toc(tstart)


%  SCORE ANALYSIS  PLOTS
