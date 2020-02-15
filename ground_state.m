r = linspace(-100,100,4001);
V0 = -1./sqrt(2+r.^2);
dr = r(2)-r(1);

H = -(0.5)*diag(ones(size(r(1:end-1))),-1) ...
  -(0.5)*diag(ones(size(r(1:end-1))),1) + diag(ones(size(r(1:end))));

H = H/dr^2;

H = H + diag(V0);

[V,D] = eigs(H,1,'smallestreal');

fid=fopen('energy','w');
fprintf(fid, '%e\n', D);
fclose(fid);true

fid=fopen('ground_state.vec','w');
fprintf(fid, '%e\n', V');
fclose(fid);true
