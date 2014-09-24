function y = resamp_lowpass(x, r)

[L, M] = rat(r);
xe = zeros(1,L*length(x));
xe(1:L:end) = x;
if (L>1 | M >1)
    h = lowpass(L, min(pi/M, pi/L));
    xi = conv(xe,h);

    xi = xi((length(h)+1)/2:end - (length(h) -1)/2);
else
    xi = xe;
end

y = xi(1:M:end);

function h=lowpass(gain,wc)

dw = 0.1*pi;
delta = 0.001;
A = -20*log10(delta);
beta = 0.1102*(A-8.7)*(A>50) + (0.5842*(A-21)^0.4 + 0.07886*(A-21))*(A<=50)*(A>=21);

N = ceil((A-8)/2.285/dw/2)*2;
h = gain*fir1(N,wc/pi,kaiser(N+1,beta));

 
