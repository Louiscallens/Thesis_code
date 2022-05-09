Nb_tries = 1000;
for m = 1:Nb_tries
problem.nu = 2;

M = mesh(10, 1, 5, 0, problem.nu);
M.Nu = NaN + zeros(problem.nu, 9);
M.Nu(:,1) = 2; M.Nu(:,2) = 1; M.Nu(:,3) = 0;
M.Nu(:,4:end) = round(2.*rand(size(M.Nu(:,4:end))));

U = cell(problem.nu, 0);
for i = 1:length(M.s)-1
    for l = 1:problem.nu
        if M.Nu(l,i) == 0 || M.Nu(l,i) == 1
            U{l}{i} = rand(1, 1);
        else
            U{l}{i} = rand(1, M.Nk(i));
        end
    end
end
for l = 1:problem.nu
    if M.Nu(l,end) ~= 0
        U{l}{end+1} = rand(1,1);
    end
end

Uex = expandU(U, M);

u = [Uex{:}];
svals = [];
for k = 1:length(M.s)-1
    svals = [svals, M.sc{k}(1:end-1)];
end
svals = [svals, M.s(end)];
if (length(svals) ~= size(u,2))
    disp('oops');
end
%figure; plot(svals, u', '.-');
%disp(M.Nu);
end