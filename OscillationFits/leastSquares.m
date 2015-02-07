function sum_square_errors = leastSquares(pp,tt,yy)
  [a,b,c,d]=pp;
  y_trial = dampedOscModel(a,b,c,d,tt);
difference = yy-y_trial;
sum_square_errors = sum(difference.^2);
