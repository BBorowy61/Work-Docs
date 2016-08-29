
avg=zeros(size(x));
W=200;
for i=1:length(x),
e=i+W-1;
if e>length(x)
    e=length(x);
end
Sum = sum(x(i:e));
avg(floor((e+i)/2))=Sum/W;
end
plot(t,x,t,avg);shg