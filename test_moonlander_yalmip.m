% test moonlander problem with a solver different from IPOPT
yalmip('clear')

%% define some method variables
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

%% define optimization variables
X = {}; U = {};
for k = 1:Nb_inter
    X{k} = sdpvar(2,5);
    U{k} = sdpvar(1,5);
end
X{end+1} = sdpvar(2,1);
tf = sdpvar(1,1);
rhs = @(x, u, t) [x(2,:); -1.5 + u];
rhsa = @(x, u, t) 2.*rhs(x, u, t)./(tf+1);

%% define constraints 
Constraints = [];

% initial and final constraint
Constraints = [Constraints; X{1}(:,1) == [10; -2]];
Constraints = [Constraints; X{end}(:,end) == [0; 0]];

% dynamics
for k = 1:Nb_inter
    tcoll = M.tc{k};
    xcoll = [X{k}, X{k+1}(:,1)];
    for j = 1:length(tcoll)-1
        diff = LagrangePolynomialDiff2(tcoll, xcoll, tcoll(j));
        Constraints = [Constraints;...
            diff == rhsa(X{k}(:,j), U{k}(:,j), tcoll(j))];
    end
end

% pathconstraints
for i = 1:Nb_inter
    Constraints = [Constraints; U{k} <= 3];
    Constraints = [Constraints; U{k} >= 0];
end

%% define an objective
Objective = tf;

%% set solver options
%options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
options = sdpsettings('verbose', 1, 'solver', 'FMINCON');

%% solve the problem
sol = optimize(Constraints,Objective,options);

%% post-processing
if sol.problem == 0
 % Extract and display value
 sol = value(x);
 xs = sol(1:(N+1)*nx);
us = sol((N+1)*nx + 1:end);

x1s = zeros(size(N+1,1));
x2s = zeros(size(N+1,1));
for i = 1:N+1
    x1s(i) = xs(2*i-1);
    x2s(i) = xs(2*i);
end
figure; hold on;
plot(x1s);
plot(x2s);
plot(us);
legend('x1', 'x2', 'u');
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end