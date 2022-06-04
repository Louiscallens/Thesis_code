classdef chicane2
    
    properties
        problem;
        R1; R2;
        s1; s2; s3; s4;
        total_length;
    end
    
    methods
        function obj = chicane2(problem)
            obj.problem = problem;
            obj.R1 = problem.L2 - problem.L1;
            obj.R2 = problem.L3 - problem.L2;

            obj.s1 = problem.L1;
            obj.s2 = obj.s1 + pi*obj.R1/2;
            obj.s3 = obj.s2 + pi*obj.R2/2;
            obj.s4 = obj.s3 + problem.L4-problem.L3;
            obj.total_length = obj.s4;
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
                elseif s <= obj.s2
                    X = obj.problem.L1 + obj.R1*sin((s-obj.s1)/(obj.s2-obj.s1)*pi/2);
                    Y = obj.R1*(1-cos((s-obj.s1)/(obj.s2-obj.s1)*pi/2));
                elseif s < obj.s3
                    X = obj.problem.L1 + obj.R1 + obj.R2*(1-cos((s-obj.s2)/(obj.s3-obj.s2)*pi/2));
                    Y = obj.R1 + obj.R2*sin((s-obj.s2)/(obj.s3-obj.s2)*pi/2);
                else
                    X = obj.problem.L1 + obj.R1 + obj.R2 + (obj.problem.L4-obj.problem.L3)*(s-obj.s3)/(obj.s4-obj.s3);
                    Y = obj.R1 + obj.R2;
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
                elseif s <= obj.s2
                    rho = obj.R1;
                elseif s < obj.s3
                    rho = -obj.R2;
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
                elseif s <= obj.s2
                    angle = pi*(s-obj.s1)/(2*(obj.s2-obj.s1));
                elseif s < obj.s3
                    angle = pi/2*(1-(s - obj.s2)/(obj.s3-obj.s2));
                else
                    angle = 0;
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
                elseif s <= obj.s2
                    der = pi/(2*(obj.s2-obj.s1));
                elseif s < obj.s3
                    der = -pi/(2*(obj.s3-obj.s2));
                else
                    der = 0;
                end
                ders(i) = der;
            end
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

