
u = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 1 ];
%u = [1 1 1 1 1 1 1 0; 0 0 0 0 0 0 0 1];
x = [1 3; 1.5 3.2; 1.3 2.8; 3 1; 4 1; 3.5 1.5];
%res = cluster_finder(u, x);
%dist = distance_finder(res, x);
%disp(dist);
%new = U_updater(dist);
%disp(new);
%disp(res);
%disp(norm_calculator(u, new));
ress = fcmean(x, 100, u, 0.01);
disp(ress);
c1 = [.5 .5 .5];
scatter(x(:,1), x(:,2), [], c1, 'filled');
c2 = [.2 .2 .2];
%disp(ress(1,:));
hold
scatter(ress(:,1), ress(:,2), [], c2);
disp(ress(:,2) + " in");
disp(ress)



function [cluster_centers] = fcmean(data_matrix, number_iterations, membership_matrix, thresh)
    while number_iterations > 0 
        loc = cluster_finder(membership_matrix, data_matrix);
        dist = distance_finder(loc, data_matrix);
        u = U_updater(dist);
        membership_matrix = u;
        number_iterations = number_iterations - 1;
    end
    cluster_centers = loc;
end


function [c_centers] = cluster_finder(initial_U, data_locations)
    [m,n] = size(initial_U);
    c_centers = zeros(m,2);
    m_param = 2;
    for i = 1:1:m
        temp1=0;
        temp2=0;
        temp3=0;
        for j = 1:1:n
            temp1 = temp1 + data_locations(j,1) * ((initial_U(i,j))^m_param);
            temp3 = temp3 + data_locations(j,2) * ((initial_U(i,j))^m_param);
            temp2 = temp2 + ((initial_U(i,j))^m_param);
        end
        c_centers(i,1) = temp1 / temp2;
        c_centers(i,2) = temp3 / temp2;
    end
end


function [distance_matrix] = distance_finder(cluster_centers, data_locations)
    [m,n] = size(cluster_centers);
    [x,y] = size(data_locations);
    distance_matrix = zeros(m, x);
    for i = 1:1:m
        cluster_x = cluster_centers(i,1);
        cluster_y = cluster_centers(i,2);
        for j = 1:1:x
           data_x = data_locations(j,1);
           data_y = data_locations(j,2);
           distance = sqrt((cluster_x - data_x)^2 + (cluster_y - data_y)^2);
           distance_matrix(i,j) = distance;
        end
    end
end


function [new_U] = U_updater(distance_centers)
    [m,n] = size(distance_centers);
    new_U = zeros(m,n);
    for i=1:1:n
        for j=1:1:m
            temp1 = distance_centers(j, i);
            temp2 = 0;
            for k=1:1:m
                temp2 = temp2 + (temp1/distance_centers(k,i))^2;
            end
            new_U(j,i) = 1 / temp2;
            if isnan(new_U(j,i))
                new_U(j,i) = 1;
            end
        end    
    end
end


function [norm_result] = norm_calculator(new_U, old_U)
    [m,n] = size(new_U);
    norm_result = 0;
    for i = 1:1:m
        for j=1:1:n
            if abs(new_U(i,j)-old_U(i,j)) > norm_result
                norm_result = abs(new_U(i,j)-old_U(i,j));
            end
        end
    end
    
end


