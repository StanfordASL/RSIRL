function [j_noise, robot_action] = sample_disturbance(disturbances, counter, p_vec)
    
if (counter == 8)
    j_noise = 4;
elseif (counter == -8)
    j_noise = 5;
else
    j_noise = find(mnrnd(1,p_vec));
end

robot_action = disturbances.list{j_noise};

end
