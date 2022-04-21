classdef smooth_hairpin
    
    properties
        total_length;
        s1; s2;
        a; b; c; m;
        A; B; C; D; E;
    end
    
    methods
        function obj = smooth_hairpin(params)
            obj.total_length = params.total_length;
            obj.s1 = params.s1; obj.s2 = params.s2;
            obj.a = params.a; obj.m = params.m;
            obj.b = 2*(obj.a-obj.m)/obj.s1;
            obj.c = obj.m - obj.b*obj.s2/2;
            obj.A = 2*log(2*obj.a)/obj.b;
            obj.B = obj.A - 2*log(abs(obj.b*obj.s1-2*obj.a))/obj.b;
            obj.C = -obj.s1/obj.m;
            obj.D = obj.B + obj.C + obj.s2/obj.m;
            obj.E = -2*log(abs(obj.b*obj.s2+2*obj.c))/obj.b;
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
                if s <= obj.s1
                    rho = obj.a - obj.b*s/2;
                elseif s <= obj.s2
                    rho = obj.m;
                else
                    rho = obj.c + obj.b*s/2;
                end
                rhos(i) = rho;
            end
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                s = svalues(i);
                if s <= obj.s1
                    angle = -2*log(abs(obj.b*s-2*obj.a))/obj.b + obj.A;
                elseif s <= obj.s2
                    angle = obj.B + obj.C + s/obj.m;
                else
                    angle = obj.D + obj.E + 2*log(abs(obj.b*s+2*obj.c))/obj.b;
                end
                angles(i) = angle;
            end
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            ders = NaN + zeros(size(svalues));
            
            for i = 1:length(svalues)
                s = svalues(i);
                if s <= obj.s1
                    der = 2/(2*obj.a-obj.b*s);
                elseif s <= obj.s2
                    der = 1/obj.m;
                else
                    der = 2/(2*obj.c+obj.b*s);
                end
                ders(i) = der;
            end
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

