pkg load optim

psiarray = load('l1/psi25.dat');
rangemin = 25;
rangemax = 300;
freq = 0.10727;
  amp = 0.2;
  phase =0.0;
  damping = -0.13;

tt = psiarray(:,1);
psi = psiarray(:,2);
%plot(tt(1:500),psi(1:500), tt, exp(-0.15.*tt))


nn=1;
for kk = 1:length(tt)
	   if ((tt(kk)>rangemin) && (tt(kk)<rangemax))
	      ttcut(nn) = tt(kk);
              psicut(nn)=psi(kk);
              nn=nn+1;
            endif
endfor
nmax=nn-1;
lnpsi=log(abs(psicut));
%plot(ttcut,lnpsi)

ytest = 7.0.*ttcut+3.0.*ttcut.^2+ttcut.^3+1.0;
p0=[0.0,8.0,2.0,1.0,1.0];
#p0=[freq,amp,phase,damping];
pp=fminsearch(@(pp) test,p0,[],[],ttcut,ytest);
[a,b,c,d,f]=pp;
#plot(tt,lnpsi,'r-',tt,dampedOscModel(a,b,c,d,tt),'g-');

