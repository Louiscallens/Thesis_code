classdef maskingMesh
   
    properties
        s       % edges of the intervals
        sc      % cell array where each cell contains the collocation points of that interval
        Nk      % List to indicate the number of collocation points in every interval
        Nu      % List to indicate the type of control in every interval (constant, linear, polynomial)
    end
    
    methods
        function obj = maskingMesh(N, total_length, Nmin, minUDegree, nu, disconts)
            if nargin == 6
                %obj.s = [0, 0.7*disconts(1), linspace(0.9*disconts(1), total_length, N-2)];
                %obj.s = [0, 50, 55, 60, 95, linspace(105, total_length, 10)];
                %obj.s = [0, 46, 52, 100, linspace(115, total_length, 7)];
                obj.s = [0, 46, 52, 100, linspace(115, total_length, 7)];
                [~,idx] = min(abs(obj.s-disconts(2)));
                obj.s(idx) = disconts(2);
            else
                %obj.s = [0, linspace(2/3*total_length, total_length, N-1)];
                %obj.s = [0, 50, 55, 60, 95, linspace(105, total_length, 10)];
                obj.s = [0, 47, 95, linspace(105, total_length, 10)];
            end
            obj.Nk = Nmin.*ones(1,length(obj.s));
            obj.Nu = minUDegree+zeros(nu,length(obj.s));
            obj.sc = obj.add_collocation_times();
        end
        
        function sc = add_collocation_times(obj)
            sc = {};
            for k = 1:length(obj.s)-1
                sc{k} = get_collocation_interval(obj.s(k), obj.s(k+1), obj.Nk(k));
            end
        end
    end
end
