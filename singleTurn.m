classdef singleTurn
    
    properties
        problem;
        R1;
        s1; s2; s3;
        total_length;
    end
    
    methods
        function obj = singleTurn(problem)
            obj.problem = problem;
            obj.R1 = problem.L2 - problem.L1;

            obj.s1 = problem.L1;
            obj.s2 = obj.s1 + pi*obj.R1/2;
            obj.s3 = obj.s2 + problem.L3-problem.L2;
            obj.total_length = obj.s3;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            %svalues = obj.total_length/2.*svalues + obj.total_length/2;
            
            Xs = NaN + zeros(size(svalues));
            Ys = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    X = obj.problem.L1*s/obj.s1;
                    Y = 0;
                elseif s < obj.s2
                    X = obj.problem.L1 + obj.R1*sin((s-obj.s1)/(obj.s2-obj.s1)*pi/2);
                    Y = obj.R1*(1-cos((s-obj.s1)/(obj.s2-obj.s1)*pi/2));
                else
                    X = obj.problem.L1 + obj.R1;
                    Y = obj.R1 + obj.problem.L3*(s-obj.s2)/(obj.s3-obj.s2);
                end
                Xs(i) = X; Ys(i) = Y;
            end
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            %svalues = obj.total_length/2.*svalues + obj.total_length/2;
            
            rhos = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    rho = 1.0e30;
                elseif s < obj.s2
                    rho = obj.R1;
                else
                    rho = 1.0e30;
                end
                rhos(i) = rho;
            end
        end
        
        function angles = evaluate_angle(obj, svalues)
            %svalues = obj.total_length/2.*svalues + obj.total_length/2;
            
            angles = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    angle = 0;
                elseif s < obj.s2
                    angle = pi*(s-obj.s1)/(2*(obj.s2-obj.s1));
                else
                    angle = pi/2;
                end
                angles(i) = angle;
            end
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            %svalues = obj.total_length/2.*svalues + obj.total_length/2;
            
            ders = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    der = 0;
                elseif s < obj.s2
                    der = pi/(2*(obj.s2-obj.s1));
                else
                    der = 0;
                end
                ders(i) = der;
            end
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(0, obj.total_length, 200);
            [Xs, Ys] = obj.evaluate_track_param(svalues);

            figure; plot(Xs, Ys); axis equal;
        end
        
        function test_evaluate_radius_curvature(obj)
            svalues = linspace(0, obj.total_length, 200);
            rhos = obj.evaluate_radius_curvature(svalues);
            figure; plot(svalues, rhos);
        end
        
        function test_evaluate_angles(obj)
            svalues = linspace(0, obj.total_length, 200);
            angles = obj.evaluate_angle(svalues);
            ders = obj.evaluate_angle_derivative(svalues);
            figure; hold on; plot(svalues, angles); plot(svalues, ders);
            legend('angle', 'ders');
        end
    end
end

