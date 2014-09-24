load stim_init_count.avg
flow = stim_init_count(4);
fhigh = stim_init_count(5);
fstep = stim_init_count(6);
name = pwd;

load stim.avg
ncorr=size(stim,1);
nt=size(stim,2);
nt2=(nt-1)/2;
t=-nt2:nt2;
minval = min(min(stim));
maxval = max(max(stim));
axisval(1)=t(1);
axisval(2)=t(nt);
axisval(3)=minval;
axisval(4)=maxval;


nb = (-1 + sqrt(1+8*ncorr))/2;


