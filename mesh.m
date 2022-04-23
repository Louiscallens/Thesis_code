classdef mesh
   
    properties
        s       % edges of the intervals
        sc      % cell array where each cell contains the collocation points of that interval
        Nk      % List to indicate the number of collocation points in every interval
        Nu      % List to indicate the type of control in every interval (constant, linear, polynomial)
    end
    
    methods
        function obj = mesh(N, total_length, Nmin, minUDegree, disconts)
            if nargin == 5
                obj.s = init_intervals_discontinuities(disconts, N, total_length);
            else
                obj.s = linspace(0, total_length, N);
            end
            obj.Nk = Nmin.*ones(N,1);
            obj.Nu = minUDegree+zeros(N,1);
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
