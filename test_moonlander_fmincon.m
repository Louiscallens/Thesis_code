% test moonlander problem with a solver different from IPOPT

%% define some method variables
global Nb_inter Nmin M rhsa
Nb_inter = 20;
Nmin = 5;

%% create mesh
Mt = linspace(-1, 1, Nb_inter+1);
M = struct('t', Mt, 'Nk', Nmin.*ones(Nb_inter,1));
tcoll = {};
for k = 1:length(M.t)-1
    tau = get_collocation_interval(M.t(k), M.t(k+1), M.Nk(k));        
    tcoll{k} = tau;
end
M.tc = tcoll;

tcolls = [];
for k = 1:Nb_inter
    tcolls = [tcolls; M.tc{k}(1:end-1)'];
end
tcolls = [tcolls; M.t(end)];

%% define system dynamics
rhs = @(v, u) [v; -1.5 + u];
rhsa = @(v, u, tf) 2.*rhs(v, u)./(tf+1);

%% define constraints 
x0 = [10; -2]; xf = [0; 0];

% pathconstraints
lb = -inf + zeros(3*(Nb_inter*Nmin+1),1);
ub =  inf + zeros(3*(Nb_inter*Nmin+1),1);
lb(end-Nb_inter*Nmin+1:end) = 0;
ub(end-Nb_inter*Nmin+1:end) = 3;

%% define an objective
Objective = @(x) x(end);

%% define initial guess
x0 = NaN + zeros(3*(Nb_inter*Nmin+1),1);
x0(1:Nb_inter*Nmin+1) = 10 - 2.*tcolls + 1/10.*tcolls.^2;
x0(Nb_inter*Nmin+2:2*(Nb_inter*Nmin+1)) = -2 + 1/5.*tcolls;
x0(203:302) = 1/5 + 1.5;
x0(end) = 10;

%% set solver options
%options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
%options = sdpsettings('verbose', 1, 'solver', 'FMINCON');

%% solve the problem
%sol = optimize(Constraints,Objective,options);
fmincon(Objective, x0, [], [], [], [], lb, ub, @(x)myCon(x));

function [c, ceq] = myCon(x)
    global Nb_inter Nmin M rhsa
    
    c = []; ceq = [];
    ceq = [ceq; x(1) - 10; x(Nb_inter*Nmin+1 + 1)-2];
    ceq = [ceq; x(Nb_inter*Nmin+1); x(2*(Nb_inter*Nmin+1))];
    for k = 1:Nb_inter
        tcoll = M.tc{k};
        xcoll = [x((k-1)*Nmin+1:k*Nmin+1), x(Nb_inter*Nmin+1 + (k-1)*Nmin+1:Nb_inter*Nmin+1 + k*Nmin+1)];
        xcoll = xcoll';
        ucoll = x(2*(Nb_inter*Nmin+1)+1:end-1);
        for j = 1:Nmin
            diff = LagrangePolynomialDiff2(tcoll, xcoll, tcoll(j));
            ceq = [ceq; diff - rhsa(xcoll(2,j), ucoll(j), x(end))];
        end
    end
    %c = [c; x(2*(Nb_inter*Nmin+1)+1:3*(Nb_inter*Nmin+1)-1) - 3];
    %c = [c; -x(2*(Nb_inter*Nmin+1)+1:3*(Nb_inter*Nmin+1)-1)];
    c = [c; -x(end)];
end