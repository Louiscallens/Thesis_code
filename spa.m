classdef spa
    
    properties
        p;
        L;
        R;
        total_length;
    end
    
    methods
        function obj = spa(p)
            L = [10, 5, 30, 20, 100, 10, 10, 10, 10, 30, 20, 10, 10, 40, 10, 10, 10, 20, 15, 10, 20, 7, 20, 10, 50, 50, 30, 20, 10];
            R = [10^30, 10, 10^30, 80, 10^30, 20, 20, 10^30, 20, 10^30, 20, 10^30, 10, 10^30, 15, 10^30, 90, 10^30, 30, 20, 10^30, 15, 10^30, 15, 10^30, 150, 40, 20, 10];
            obj.p = 1;
            obj.L = L;
            obj.R = R;
            obj.total_length = sum(L);%p.total_length;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            x = @(s) integral(@(t) cos(obj.evaluate_angle(t)), 0, s, 'arrayvalued', true);
            y = @(s) integral(@(t) sin(obj.evaluate_angle(t)), 0, s, 'arrayvalued', true);
            
            Xs = NaN + zeros(size(svalues)); 
            Ys = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                Xs(i) = x(svalues(i));
                Ys(i) = y(svalues(i));
            end
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)           
            rhos = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                idx = 1;
                L_cumm = 0;
                while idx <= length(obj.L) && s >= L_cumm + obj.L(idx)
                    idx = idx + 1;
                    L_cumm = L_cumm + obj.L(idx-1);
                end
                rhos(i) = obj.R(min(idx, length(obj.R)));
            end
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = NaN + zeros(size(svalues));
            f = @(s) integral(@(t) obj.evaluate_angle_derivative(t), 0, s, 'arrayvalued', true);
            for i = 1:length(svalues)
                angles(i) = f(svalues(i));
            end
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            rhos = obj.evaluate_radius_curvature(svalues);
            ders = 1./rhos;
            ders(ders < 1e-25) = 0;
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(0, obj.total_length, 10);
            [Xs, Ys] = obj.evaluate_track_param(svalues);

            figure; plot(Xs, Ys); axis equal;
        end
        
        function test_evaluate_radius_curvature(obj)
            svalues = linspace(0, obj.total_length, 1000);
            rhos = obj.evaluate_radius_curvature(svalues);
            figure; plot(svalues, rhos);
        end
        
        function test_evaluate_angles(obj)
            svalues = linspace(0, obj.total_length, 1000);
            %svalues = [obj.L(1) + 0.01, obj.L(1) + 1];
            angles = obj.evaluate_angle(svalues);
            ders = obj.evaluate_angle_derivative(svalues);
            figure; hold on; plot(svalues, angles); plot(svalues, ders);
            legend('angle', 'ders');
        end
    end
end

