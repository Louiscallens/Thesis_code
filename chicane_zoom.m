classdef chicane_zoom
    
    properties
        problem;
        total_length;
        s_start; s_end;
        complete_track;
    end
    
    methods
        function obj = chicane_zoom(s_start, s_end, problem)
            obj.problem = problem;
            obj.s_start = s_start; obj.s_end = s_end;
            obj.complete_track = chicane(problem);
            obj.total_length = obj.s_end - obj.s_start;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            [Xs, Ys] = obj.complete_track.evaluate_track_param(obj.s_start + svalues);
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = obj.complete_track.evaluate_radius_curvature(obj.s_start + svalues);
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = obj.complete_track.evaluate_angle(obj.s_start + svalues);
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = obj.complete_track.evaluate_angle_derivative(obj.s_start + svalues);
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

