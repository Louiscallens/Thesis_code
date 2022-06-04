load('workspaces/thesis/final_problem3/high_accuracy_hp/iteration_5.mat');
string = "";
for i = 1:length(timings)
    string = string + num2str(i) + " & " + num2str(timings(i)) + " & " + iterCounts(i) + "\\\\ \n";
end
string = string + " \\hline Total & " + num2str(sum(timings)) + " & " + num2str(sum(iterCounts)) + " \\\\ \\hline";
disp(sprintf(string));