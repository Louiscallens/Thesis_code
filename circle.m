classdef circle
    
    properties
        total_length;
        circ;
        R;
    end
    
    methods
        function obj = circle(params)
            obj.total_length = params.total_length;
            obj.R = params.R;
            obj.circ = 2*pi*obj.R;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            Xs = obj.R.*sin(svalues./obj.circ*2*pi);
            Ys = obj.R.*(1-cos(svalues./obj.circ*2*pi));
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = obj.R.*ones(size(svalues));
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = svalues./obj.circ.*2*pi;
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = 2*pi/obj.circ.*ones(size(svalues));
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

