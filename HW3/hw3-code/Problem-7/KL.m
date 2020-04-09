function KL

import brml.*

% Load p.mat
load p
p=setpotclass(p,'array');
X=3; Y=3; H=3;
[x y h]=assign(1:3);

qxy = array([x y],normp(rand([X Y])));
qh = array(h,normp(rand(H,1)));

% using a general mean field update approach:
for loop=1:100
    qxy = condpot(exppot(sumpot(multpots([logpot(p) qh]),h)));
    qh = condpot(exppot(sumpot(multpots([logpot(p) qxy]),[x y])));
    qxyh = multpots([qxy qh]);
    kl(loop) = KLdiv(qxyh,p);
    plot(kl); title('KL divergence'); drawnow
end

approx = table(multpots([qxy qh]));
figure

% Print the result
fprintf(1,'The value of the minimal KL Divergence for the optimal q = %g\n', mean(abs(approx(:)-p.table(:))));