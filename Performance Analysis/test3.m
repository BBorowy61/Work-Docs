CorrelSig2Respsec = zeros(H1int, M);
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