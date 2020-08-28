
stepsize_check_convergence = ts;
increments = floor(ts/stepsize_check_convergence); 
% Now get time dependent rate matrix in and out.
% Tag when transitions occur  
num_occupied_states = length(states_num);
states_shifted = [states_in_ts(2:ts);states_in_ts(1)];
difference = states_shifted - states_in_ts;
trans_prob = zeros(num_occupied_states);
trans_prob_bytime = zeros(num_occupied_states, num_occupied_states*increments);
trans_rate_bytime = zeros(num_occupied_states, num_occupied_states*increments);
equil_dist = 0;
index_ff = zeros(1,increments);
for m = 1:increments
    n = m*stepsize_check_convergence;
    j = n - (stepsize_check_convergence - 1);

    trans_prob = trans_prob.*0;
    % Count transitions from i->j
    for i = j:n-1        
        oldstate = states_in_ts(i);
        newstate = difference(i) + states_in_ts(i);
        trans_prob(oldstate,newstate) = trans_prob(oldstate,newstate) + 1;
    end
    % Divide by time spent in i
    time_in_oldstate = sum(trans_prob,2);
    num_trans = sum(sum(trans_prob));
    trans_rate = trans_prob./num_trans;
    for i = 1:length(time_in_oldstate)
        if time_in_oldstate(i)>0
           trans_prob(i,:) = trans_prob(i,:)./time_in_oldstate(i);
        end
    end
    for i = 1:num_occupied_states
        trans_rate(i,i) = 0;
        trans_rate(i,i) = -1*sum(trans_rate(i,:));
    end   
    p = m*num_occupied_states - num_occupied_states + 1;
    q = m*num_occupied_states;
    trans_prob_bytime(:,p:q) = trans_prob; 
    trans_rate_bytime(:,p:q) = trans_rate;
    
    I = eye(num_occupied_states);
    P_minus_I = trans_prob - I;

    lhs = [P_minus_I'; ones(1,num_occupied_states)];
    rhs = [zeros(num_occupied_states, 1); 1]';
    rhs = rhs';

    equil_dist_int = lhs\rhs;
    equil_dist = equil_dist + equil_dist_int;
end
equil_dist = equil_dist./increments;
equil_dist_tomult = equil_dist';

% This is done twice only to account for the total trajectory lenght
% not being a multiple of increment size. It will not matter if icrement
% size is "ts".

trans_prob = trans_prob.*0;
    % Count transitions from i->j
for i = 2:ts        
        oldstate = states_in_ts(i);
        newstate = difference(i) + states_in_ts(i);
        trans_prob(oldstate,newstate) = trans_prob(oldstate,newstate) + 1;
end
    % Divide by time spent in i
time_in_oldstate = sum(trans_prob,2);
num_trans = sum(sum(trans_prob));
trans_rate = trans_prob./num_trans; 
for i = 1:length(time_in_oldstate)
        if time_in_oldstate(i)>0
           trans_prob(i,:) = trans_prob(i,:)./time_in_oldstate(i);
        end
end
for i = 1:num_occupied_states
        trans_rate(i,i) = 0;
        trans_rate(i,i) = -1*sum(trans_rate(i,:));
end   

I = eye(num_occupied_states);
P_minus_I = trans_prob - I;

 lhs = [P_minus_I'; ones(1,num_occupied_states)];
 rhs = [zeros(num_occupied_states, 1); 1]';
 rhs = rhs';

equil_dist_total = lhs\rhs;
    
 figure
 plot(state_prob(:,2), '-ob')
 hold on
 plot(equil_dist_total, '-xr')
 
[rel_prob, eigenv] = eig(trans_prob');
[vec] = real(eigenv)>.999;

[ind1,ind2] = find(vec==1);
equil_dist_eig = rel_prob(:,ind1)/sum(rel_prob(:,ind1));

