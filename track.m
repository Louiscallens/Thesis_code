classdef track
    
    properties
        total_length;
    end
    
    methods
        function obj = track(params)
            obj.total_length = 0;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            Xs = 0; Ys = 0;
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = 0;
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = 0;
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = 0;
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(0, obj.total_length, 1000);
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

