% Average 1st shell occupancy as function of minimum Voronoi border l
at2 = readtable('rts_fs.txt');
at2_dppc_fs = at2{:,1};
ts = length(at2_dppc_fs)
time = 1:1:ts;
time = 15000*time./max(time); % adapt for time between frames; "time" is just a vector indicating time at each frame.
time = time./1000;
at2_chol_fs = at2{:,2};
at2_dopc_fs = at2{:,3};
at2_dppc_ss = at2{:,4};
at2_chol_ss = at2{:,5};
at2_dopc_ss = at2{:,6};
at2_fs = at2_dppc_fs + at2_chol_fs + at2_dopc_fs;
at2_fs = mean(at2_fs);
at2_fs_whole = [at2_dppc_fs, at2_chol_fs, at2_dopc_fs];
at2pt_var_dppc = var(at2_dppc_fs);
at2pt_var_chol = var(at2_chol_fs);
at2pt_var_dopc = var(at2_dopc_fs);

% Compile statistics, only for border length 2 Ang.
fs_occupancy_by_lip = [at2_dppc_fs, at2_chol_fs,at2_dopc_fs ];
fs_occupancy_by_lip_avg_2Ang = sum(fs_occupancy_by_lip,1)/ts
fs_occupancy_avg_tot_2Ang = sum(fs_occupancy_by_lip_avg_2Ang,2)
fs_occupancy = sum(fs_occupancy_by_lip,2);

% Time series of 1st shell occupancy
% First compute total 1st shell occupancy for each edge length
at2_fs = at2_dppc_fs + at2_chol_fs + at2_dopc_fs;

% Rounding scheme and possible occupied states
percent_dppc_int = ones(10000,1);
for i = 1:1:100,
    q = (i - 1)*100 + 1;
    percent_dppc_int(q:q+99) = ones(100,1)*i;
end
percent_dppc = ones(1000000,1);
for i = 1:100,
    q = (i - 1)*10000 + 1;
    percent_dppc(q:q+9999) = percent_dppc_int;
end
percent_chol_int = 1:100;
percent_chol_int = percent_chol_int';
percent_chol = ones(1000000,1);
for i = 1:10000,
    q = (i - 1)*100 + 1;
    percent_chol(q:q+99) = percent_chol_int;
end
percent_dopc = ones(1000000,1);
for i = 1:100,
    q = (i - 1)*10000 + 1;
    percent_dopc(q:q+9999) = ones(10000,1)*i;
end
tot = percent_dppc + percent_chol + percent_dopc;
Tot = tot;
Tot(tot>100) = 0;
Tot(tot<100) = 0;
[fs_comp_ele] = find(Tot==100);
% fs_comp are the bins of 1st shell percentage.
fs_comp = [percent_dppc(fs_comp_ele), percent_chol(fs_comp_ele), percent_dopc(fs_comp_ele)];

% Now compute 1st shell percentage from at2_fs_whole and at2_fs:
at2_fs = at2_dppc_fs + at2_chol_fs + at2_dopc_fs;
percent_dppc_column = at2_fs_whole(:,1)./at2_fs*100;
percent_chol_column = at2_fs_whole(:,2)./at2_fs*100;
percent_dopc_column = at2_fs_whole(:,3)./at2_fs*100;

% Change rounding scheme here. The correction for number of states will be
% automatic.
percent_dppc_column_int = round(percent_dppc_column./5)*5;
percent_chol_column_int = round(percent_chol_column./5)*5;
percent_dopc_column_int = round(percent_dopc_column./5)*5;

sumofthem = percent_dppc_column_int + percent_chol_column_int + percent_dopc_column_int;

all_column = [percent_dppc_column, percent_chol_column, percent_dopc_column];
for i = 1:ts,
    if sumofthem(i) > 100, 
      a = 10*rem(percent_dppc_column(i)*.1,1); 
      a = round(a);
      b = 10*rem(percent_chol_column(i)*.1,1);
      b = round(b);
      c = 10*rem(percent_dopc_column(i)*.1,1);
      c = round(c);
      d = [a,b,c];
      % Say if it equals 3 or 8. Then again if it equals 8 or 9. 
      ind = find(d==3|d==8|d==4|d==9);
      ind2 = find(min(d(ind)));
      if length(ind) == 2,
          if ind ==[2 3], 
          ind2 = ind2+1;
          end
          if ind == [1 3];
             if ind2 == 2,
               ind2 = 3;
             end
          end
      end
      all_column(i,ind2) = floor(all_column(i,ind2)./5)*5;
    end
end

for i = 1:ts,
    if sumofthem(i) < 100, 
      a = 10*rem(percent_dppc_column(i)*.1,1); 
      a = round(a);
      b = 10*rem(percent_chol_column(i)*.1,1);
      b = round(b);
      c = 10*rem(percent_dopc_column(i)*.1,1);
      c = round(c);
      d = [a,b,c];
      ind = find(d==1|d==2|d==6|d==7);
      ind2 = find(max(d(ind)));
      if length(ind) == 2,
          if ind ==[2 3], 
          ind2 = ind2+1;
          end
          if ind == [1 3];
             if ind2 == 2,
               ind2 = 3;
             end
          end
      end
      all_column(i,ind2) = ceil(all_column(i,ind2)./5)*5;
    end
end
      

per_dppc_column = round(all_column(:,1)./5)*5;
per_chol_column = round(all_column(:,2)./5)*5;
per_dopc_column = round(all_column(:,3)./5)*5;


percent_fs_comp = [per_dppc_column,per_chol_column,per_dopc_column];
unique_states = unique( percent_fs_comp(:,[1 2 3]), 'rows', 'stable');

[filled_fs_comp, indA, indB] = intersect(fs_comp, unique_states, 'rows', 'stable'); 

% Label states with number and find state probabilities by number.
% Create matrix of state labels.
states_num = 1:1:length(filled_fs_comp(:,1));
states_num = states_num';
states_label = [filled_fs_comp, states_num];

% Flag state at each timestep in traj and compute state probabilities each
% *stepsize_check_convergence steps.
stepsize_check_convergence = ts; % You can make this smaller to break your trajectory into chunks for statistical analylsis.
increments = floor(ts/stepsize_check_convergence);
states_in_ts = ones(ts,1);
freq = zeros(length(states_num),increments);
state_prob = freq;
for i = 1:length(states_num),
    for j = 1:increments,
        q = j*stepsize_check_convergence;
        [ind3,~] = ismember(percent_fs_comp(1:q,:),filled_fs_comp(i,:),'rows');
        freq(i,j) = sum(ind3);
    end
    j = ts;
    [ind3,ind4] = ismember(percent_fs_comp(1:j,:),filled_fs_comp(i,:),'rows');      
    ind3 = ind3*i; 
    states_in_ts_int = states_in_ts.*ind3;
    states_in_ts = states_in_ts + states_in_ts_int;
end

% This is the list of states by timestep
states_in_ts = states_in_ts - 1;

% Matrix of state probabilities of form (state#, prob) for each increment
% in time.
for i = 1:increments,
    q = i*stepsize_check_convergence;
    state_prob(:,i) = freq(:,i)/q;
end
    state_prob = [states_num, state_prob];