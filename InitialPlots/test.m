function pf = test(pp,tt,yy)
  [a,b,c,d,e]=pp;

y_trial = a.*tt+b.*tt.^2+c.*tt.^3*d;
difference = yy-y_trial;
testO = sum(difference.^2);
