classdef final_problem
    
    properties
        L; Ls;
        R; rx; ry;
        angles;
        Xs; Ys;
        total_length;
        a_init;
    end
    
    methods
        function obj = final_problem()
            obj.R =    [1e30, 20,    15, 1e30,   5, 1e30,   5, 1e30, 60, 1e30,  18, 1e30];
            L =        [100,  NaN,  NaN,   50, NaN,    5, NaN,   107.87, NaN,   5.75, NaN, 30];
            rotation = [NaN,  250, -250,  NaN,  90,  NaN,  90,  NaN,  60,  NaN, 120, NaN];
            obj.a_init = 0;
            %obj.R =    [1e30, 5,  1e30, 300, 1e30,  10,  10, 1e30,   10, 1e30,   10, 1e30,  15, 1e30,   20, 1e30,   40, 1e30,   10,  20, 1e30,   15, 1e30,   20, 1e30, 50,  10,  10, 1e30];
            %L =        [20,   NaN,  30, NaN,  100, NaN, NaN,   10,  NaN,   20,  NaN,   10, NaN,   30,  NaN,   10,  NaN,   30,  NaN, NaN,   20,  NaN,   15,  NaN,   78, NaN, NaN, NaN, 5];
            %rotation = [NaN, -160, NaN, -10,  NaN, -90,  90,  NaN,  -90,  NaN, -180,  NaN,  90,  NaN,   90,  NaN,   45,  NaN, -100,  70,  NaN,  -90,  NaN,  -100,  NaN,  100, -105,  80, NaN];
            %obj.a_init = pi;
            [obj.R, obj.L] = obj.merge_lengths_rotation(L, rotation);
            obj.Ls = obj.compute_cummulated_L(); 
            obj.angles = obj.compute_angles();
            [obj.rx, obj.ry, obj.Xs, obj.Ys] = obj.compute_coordinates();
            obj.total_length = sum(obj.L);%p.total_length;
        end
        
        function [R, L] = merge_lengths_rotation(obj, L, rotation)
            R = obj.R;
            for i = 1:length(L)
                if R(i) < 1e25
                    L(i) = abs(2*pi*obj.R(i)*rotation(i)/360);
                    if rotation(i) < 0
                        R(i) = -R(i);
                    end
                end
            end
        end
        
        function Ls = compute_cummulated_L(obj)
            Ls = zeros(size(obj.L));
            Ls(1) = obj.L(1);
            sum = obj.L(1);
            for i = 2:length(Ls)
                sum = sum + obj.L(i);
                Ls(i) = sum;
            end
        end
        
        function angles = compute_angles(obj)
            angles = zeros(size(obj.L));
            angles(1) = obj.a_init;
            for i = 2:length(angles)
                if obj.R(i) > 1e25
                    angles(i) = angles(i-1);
                else
                    angles(i) = angles(i-1) + obj.L(i)/obj.R(i);
                end
            end
        end
        
        function [rx, ry, Xs, Ys] = compute_coordinates(obj)
            Xs = zeros(size(obj.L)); Ys = zeros(size(obj.L));
            rx = NaN + zeros(size(obj.L)); ry = NaN + zeros(size(obj.L));
            Xs(1) = cos(obj.angles(1))*obj.L(1);
            Ys(1) = sin(obj.angles(1))*obj.L(1);
            for i = 2:length(Xs)
                if obj.R(i) > 1e25
                    Xs(i) = Xs(i-1) + cos(obj.angles(i-1))*obj.L(i);
                    Ys(i) = Ys(i-1) + sin(obj.angles(i-1))*obj.L(i);
                else
                    rx(i) = Xs(i-1) - sin(obj.angles(i-1))*obj.R(i);
                    ry(i) = Ys(i-1) + cos(obj.angles(i-1))*obj.R(i);
                    
                    phi = obj.L(i)/obj.R(i);
                    rot = [cos(phi), -sin(phi); sin(phi), cos(phi)];
                    new = rot*[Xs(i-1)-rx(i); Ys(i-1)-ry(i)] + [rx(i); ry(i)];
                    Xs(i) = new(1);
                    Ys(i) = new(2);
                end
            end
        end
        
        function [Xs, Ys] = evaluate_track_param(obj, svalues)            
            Xs = NaN + zeros(size(svalues));
            Ys = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                idx = 1;
                while idx < length(obj.Ls) && s >= obj.Ls(idx)
                    idx = idx + 1;
                end
                idx = min(idx, length(obj.Ls));
                
                x1 = 0 + (idx > 1)*obj.Xs(max(1,idx-1)); y1 = 0 + (idx > 1)*obj.Ys(max(1,idx-1));
                a1 = obj.angles(1) + (idx > 1)*(obj.angles(max(1,idx-1))-obj.angles(1));
                Ls1 = 0 + (idx > 1)*obj.Ls(max(1,idx-1));
                
                if obj.R(idx) > 1e25
                    Xs(i) = x1 + (s-Ls1).*cos(a1);
                    Ys(i) = y1 + (s-Ls1).*sin(a1);
                else
                    Xs(i) = obj.rx(idx) + obj.R(idx)*sin(a1+(s-Ls1)/obj.R(idx));
                    Ys(i) = obj.ry(idx) - obj.R(idx)*cos(a1+(s-Ls1)/obj.R(idx));
                end
            end
        end
        
        function rhos = evaluate_radius_curvature(obj, svalues)           
            rhos = NaN + zeros(size(svalues));

            for i = 1:length(svalues)
                s = svalues(i);
                idx = 1;
                L_cumm = 0;
                while idx <= length(obj.L) && s >= L_cumm + obj.L(idx)
                    idx = idx + 1;
                    L_cumm = L_cumm + obj.L(idx-1);
                end
                rhos(i) = obj.R(min(idx, length(obj.R)));
            end
        end
        
        function angles = evaluate_angle(obj, svalues)
            angles = interp1([0,obj.Ls], [obj.a_init,obj.angles], svalues);
        end
        
        function ders = evaluate_angle_derivative(obj, svalues)
            rhos = obj.evaluate_radius_curvature(svalues);
            ders = 1./rhos;
            ders(ders < 1e-25) = 0;
        end
        
        function test_evaluate_track_param(obj)
            svalues = linspace(0, obj.total_length, 1000);
            [Xs, Ys] = obj.evaluate_track_param(svalues);

            %figure;
            clf; plot(Xs, Ys); axis equal; hold on;
            plot(obj.Xs, obj.Ys, '.b');
        end
        
        function test_evaluate_radius_curvature(obj)
            svalues = linspace(0, obj.total_length, 100);
            rhos = obj.evaluate_radius_curvature(svalues);
            figure; plot(svalues, rhos);
        end
        
        function test_evaluate_angles(obj)
            svalues = linspace(0, obj.total_length, 100);
            angles = obj.evaluate_angle(svalues);
            ders = obj.evaluate_angle_derivative(svalues);
            figure; hold on; plot(svalues, angles); plot(svalues, ders);
            legend('angle', 'ders');
        end
    end
end

