classdef maskingMesh2
   
    properties
        s       % edges of the intervals
        sc      % cell array where each cell contains the collocation points of that interval
        Nk      % List to indicate the number of collocation points in every interval
        Nu      % List to indicate the type of control in every interval (constant, linear, polynomial)
    end
    
    methods
        function obj = maskingMesh2(N, total_length, Nmin, minUDegree, nu, disconts)
            obj.s = [20*linspace(0,1,4), 25:15:100];
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
