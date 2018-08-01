function paths = create_2D_path(n_paths)

%plunk down some points in 2D space
n_nodes = 5;

paths = cell(n_paths,1);
for i = 1:n_paths
    path = zeros(n_nodes,2);
    path(1,:) = [-3, -3] + [6,6].*rand(1,2);
    for j = 2:n_nodes
            path(j,:) = [-2, -2] + [4, 4].*rand(1,2) + path(j-1,:);
    end
    paths{i} = path;
end

figure()
title('waypoint path');
hold on
for k = 1:n_paths
    for j = 1:n_nodes-1
        line([paths{k}(j,1), paths{k}(j+1,1)],[paths{k}(j,2),paths{k}(j+1,2)],...
            'Color', 'b', 'LineWidth', 2);
    end
end

end
