classdef generic_track
    
    properties
        total_length;
        rho;
    end
    
    methods
        function obj = generic_track(params)
            obj.rho = params.rho;
            obj.total_length = params.total_length;
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
            rhos = obj.rho(svalues);
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = NaN + zeros(size(svalues));
            f = @(s) integral(@(t) obj.evaluate_angle_derivative(t), 0, s, 'arrayvalued', true);
            for i = 1:length(svalues)
                angles(i) = f(svalues(i));
            end
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = 1./obj.rho(svalues);
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(0, obj.total_length, 100);
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
            angles = obj.evaluate_angle(svalues);
            ders = obj.evaluate_angle_derivative(svalues);
            figure; hold on; plot(svalues, angles); plot(svalues, ders);
            legend('angle', 'ders');
        end
    end
end

