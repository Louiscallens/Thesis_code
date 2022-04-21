classdef smooth_track
    
    properties
        total_length;
        a;
        b;
    end
    
    methods
        function obj = smooth_track(problem)
            obj.total_length = 100;
            obj.a = 10;
            obj.b = 0.1;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            Xs = svalues;
            Ys = -obj.a+obj.a.*cos(obj.b.*svalues);
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = -min(max(1/(obj.a*obj.b^2)./cos(obj.b.*svalues), -1.0e30), 1.0e30);
            %rhos = -min(max(10./cos(0.1.*svalues), -1.0e30), 1.0e30);
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = -(obj.a*obj.b).*sin(obj.b.*svalues);
            %angles = -sin(0.1.*svalues);
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = -(obj.a*obj.b*obj.b).*cos(obj.b.*svalues);
            %ders = -0.1.*cos(0.1.*svalues);
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