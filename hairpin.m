classdef hairpin
    
    properties
        total_length;
        s1; s2;
        L1; R; L2;
    end
    
    methods
        function obj = hairpin(params)
            obj.L1 = params.L1;  % straight line before start of hairpin
            obj.R = params.R;   % radius of curvature of the hairpin
            obj.L2 = params.L2; % straight line after hairpin
            obj.total_length = obj.L1 + pi*obj.R + obj.L2;
            obj.s1 = obj.L1;
            obj.s2 = obj.s1 + pi*obj.R;
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)
            Xs = NaN + zeros(size(svalues));
            Ys = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    X = s;
                    Y = 0;
                elseif s < obj.s2
                    X = obj.s1 + obj.R*sin((s-obj.s1)/(obj.s2-obj.s1)*pi);
                    Y = obj.R - obj.R*cos((s-obj.s1)/(obj.s2-obj.s1)*pi);
                else
                    X = obj.s1 - (s-obj.s2);
                    Y = 2*obj.R;
                end
                Xs(i) = X; Ys(i) = Y;
            end
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)
            rhos = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    rho = 1.0e30; 
                elseif s < obj.s2
                    rho = obj.R;
                else
                    rho = 1.0e30;
                end
                rhos(i) = rho;
            end
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    angle = 0;
                elseif s < obj.s2
                    angle = pi*(s-obj.s1)/(obj.s2-obj.s1);
                else
                    angle = pi;
                end
                angles(i) = angle;
            end
                
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                s = svalues(i);
                if s < obj.s1
                    der = 0;
                elseif s < obj.s2
                    der = pi/(obj.s2-obj.s1);
                else
                    der = 0;
                end
                ders(i) = der;
            end
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

