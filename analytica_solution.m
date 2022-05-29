L = 100; amax = 20; amin = -5; v0 = 10;
c = v0^2/2; d = v0^2/2-L*amin;
s1 = -L*amin/(amax-amin);

u = @(s) (s < s1).*amax + (s >= s1).*amin;
v = @(s) (s < s1).*sqrt(2*(c+amax.*s)) + (s >= s1).*sqrt(2*(d+s*amin));

svals = linspace(0, L, 1000);

t = @(s) integral(@(s) 1./v(s), 0, s);
tvals = zeros(size(svals));
for i = 1:length(svals)
    tvals(i) = t(svals(i));
end

figure;
plot(tvals, v(svals));

figure;
plot(svals, v(svals));