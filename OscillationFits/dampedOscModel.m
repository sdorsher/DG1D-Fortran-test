function dampedOsc = dampedOscModel(freq,amp,phase,damping,time)
% damped harmonic oscillator in abs log scale
pi=3.1415926
  sho=amp*cos(2.0.*pi.*freq.*time+phase);
decay=exp(-1.0.*damping.*time);
dampedOsc=log(abs(sho.*decay));

endfunction






