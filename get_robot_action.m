function robot_action = get_robot_action(disturbances, j, K, counter)

last = disturbances.last;
if disturbances.last == 4 % robot finishing turn
    if counter == -1 % robot finishing left turn
        last = 4;
    else % robot finishing right turn
        last = 5;
    end
end

if (j == 4) % if robot turns
    if (counter == -1) % robot on left lane
        j = 5; % robot turns right
    end
end
robot_action_new = disturbances.list{j};
robot_action = [disturbances.list{last}(:,disturbances.overlap+1:K),...
    robot_action_new(:,1:disturbances.overlap)];
end
