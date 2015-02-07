psiarray= load('l2/psi25.dat');
tt = psiarray(:,1);
psi = psiarray(:,2);

logtt=log(tt);
for kk=1:length(tt)-1
	 difft(kk)=tt(kk+1)-tt(kk);
difflogt(kk) = logtt(kk+1)-logtt(kk);
endfor

title("l=2,order=25");
xlabel("Step");
ylabel("Difference in t");
plot(kk,difft,'-')
figure 
title("l=2,order=25");
xlabel("Step");
ylabel("Difference in log(t)");
plot(kk,difflogt,'-')
