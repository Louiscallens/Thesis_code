classdef straight
    
    properties
        problem;
        R1; R2;
        s1; s2; s3; s4;
        total_length;
    end
    
    methods
        function obj = straight(problem)
            obj.problem = problem;
            obj.total_length = problem.total_length;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            Xs = svalues;
            Ys = zeros(size(svalues));
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = 1.0e30 + zeros(size(svalues));
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = zeros(size(svalues));
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = zeros(size(svalues));
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(-1, 1, 100);
            [Xs, Ys] = obj.evaluate_track_param(svalues);

            figure; plot(Xs, Ys); axis equal;
        end
        
        function test_evaluate_radius_curvature(obj)
            svalues = linspace(-1, 1, 100);
            rhos = obj.evaluate_radius_curvature(svalues);
            figure; plot(svalues, rhos);
        end
        
        function test_evaluate_angles(obj)
            svalues = linspace(1, -1, 1000);
            angles = obj.evaluate_angle(svalues);
            ders = obj.evaluate_angle_derivative(svalues);
            figure; hold on; plot(svalues, angles); plot(svalues, ders);
            legend('angle', 'ders');
        end
    end
end

