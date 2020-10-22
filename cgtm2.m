% cgtm2.m populates a transition matrix "trans_prob" from states_in_ts and solves the matrix equation
% to get the equilibrium distribution of states "equil_dist_total."
% The matrix equation is: equil_dist_total = equil_dist_total*trans_prob

% The code has 3 parts.
%    1. Given total number of states "num_occupied_states", create a matrix of size num_occpuied_states x num_occpupied_states.
         % rows of this matrix represent transitions from state i to state j. For example, the first row represents transtions from
         % state 1 to all other states, and element (1,6) of the matrix corresponds to transtions from state 1 to state 6. 
%        % Populate this matrix with transition counts. 
%    2. Normalize by total time spent in state i; that is, divide each row by the sum of that row. Each row sums to 1.
%    3. Solve the matrix equation: equil_dist_total = equil_dist_total*trans_prob
%          The code does this in two ways and compares.

stepsize_check_convergence = ts;  % In this version of the code, the entire trajectory is used. 
increments = floor(ts/stepsize_check_convergence);   
num_occupied_states = length(states_num);
trans_prob = zeros(num_occupied_states);  % Initializes a variable for the transition matrix.
trans_prob_bytime = zeros(num_occupied_states, num_occupied_states*increments);
trans_rate_bytime = zeros(num_occupied_states, num_occupied_states*increments);
equil_dist = 0;
index_ff = zeros(1,increments); 
for m = 1:increments
    n = m*stepsize_check_convergence;
    j = n - (stepsize_check_convergence - 1);

    trans_prob = trans_prob.*0; % set all transition probabilities to zero before populating matrix
    % Count transitions from i->j
    for i = j:n-1        
        oldstate = states_in_ts(i);
        newstate = states_in_ts(i+1) + states_in_ts(i);
        trans_prob(oldstate,newstate) = trans_prob(oldstate,newstate) + 1;
    end
    % Divide by time spent in i
    time_in_oldstate = sum(trans_prob,2);  % This sums along the rows. Produces a vector of row sums.
    num_trans = sum(sum(trans_prob));  % This is and "trans_rate" are used later to get a time-dependent rate matrix. If you
    % sample at regular intervals, these methods give the same results.
    trans_rate = trans_prob./num_trans;
    
    % Below normalizes the rows of the transition matrix. 
    for i = 1:length(time_in_oldstate)
        if time_in_oldstate(i)>0
           trans_prob(i,:) = trans_prob(i,:)./time_in_oldstate(i);
        end
    end
    
    % Below sets the diagonals of the time depedent rate matrix to -1*sum(other elements in the row) to normalize rows of the 
    % time depednent rate matrix
    for i = 1:num_occupied_states
        trans_rate(i,i) = 0;
        trans_rate(i,i) = -1*sum(trans_rate(i,:));
    end   
    
    % I and P_minus_I are needed to create a system of equations to solve the matrix equation: 
        % equil_dist_total = equil_dist_total*trans_prob
        % I is a matrix the same size as trans_prob, with the diagonals equal to 1.
    I = eye(num_occupied_states);
    P_minus_I = trans_prob - I;
        % An additional row is added to P_minus_I to ensure the probabilities sum to 1 (they are normalized).
    lhs = [P_minus_I'; ones(1,num_occupied_states)];
    rhs = [zeros(num_occupied_states, 1); 1]';
    rhs = rhs';

    equil_dist = lhs\rhs;
end
equil_dist_tomult = equil_dist'; % This is so you can check that your result works by multiplying: equil_dist*trans_prob to
   % reproduce the equil_dist

equil_dist_total = lhs\rhs;
    
 figure
 plot(state_prob(:,2), '-ob')
 hold on
 plot(equil_dist_total, '-xr')
 
 % This is an alternate way to solve the matrix equation: 
   % The eigenvalues of the transition matrix = relative state probabilities. 
   % So if these are normalized, you should get the equil_dist.
   
[rel_prob, eigenv] = eig(trans_prob');
[vec] = real(eigenv)>.999;

[ind1,ind2] = find(vec==1);
equil_dist_eig = rel_prob(:,ind1)/sum(rel_prob(:,ind1));
% equil_dist_eig should be equil to equil_dist_tot.

