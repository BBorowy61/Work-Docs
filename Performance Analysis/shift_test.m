sh=2;
for i=1:length(aa),
    if (i+sh)>length(aa),
        aa1(i) = 0;
    else aa1(i)=aa(i+sh);
    end
end