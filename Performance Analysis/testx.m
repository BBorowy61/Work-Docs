% [HE] PrecisScoreAREG
PrecisScoreAREG = zeros(H1int,1); 
temp = Vshift(UnitResp,1);
for i = 1:H1int,
    if FleetTREG(i) ~= 0,
        PrecisScoreAREG(i) = min(max(1-abs(temp(i) - RawUnitReg(i))/UnitAREG(i),0),1);
    end
end

temp1 = Vshift(RawUnitReg,i-1);