tcoll = [0, 1, 2, 3];
xcoll = [1, 1, 1, 1; 1, 2, 3, 4; 8, 7, 5, 5];

f =  @(t) LagrangePolynomialEval(tcoll, xcoll, t);
df = @(t) LagrangePolynomialDiff2(tcoll, xcoll, t);

tvalues  = linspace(tcoll(1)-1, tcoll(end)+1, 100);
fvalues  = f(tvalues);
dfvalues = df(tvalues);

figure;
plot(tvalues, fvalues); hold on;
plot(tvalues, dfvalues);
legend('x1', 'x2', 'x3', 'dx1', 'dx2', 'dx3');