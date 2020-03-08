import brml.*

%% Initialize the setting
D=20;S=40; dv=1:D; sv=D+1:D+S;
diseaseNet = load('diseaseNet.mat')
pot=str2cell(setpotclass(pot,'array'));

%% Junction tree algorithm
tstart=tic;

% Build junction tree
[jtpot jtsep infostruct]=jtree(pot);

% Absorption
jtpot=absorption(jtpot,jtsep,infostruct);

t=toc(tstart);


for i=1:length(jtpot)
    label{i}=horzcat(variable(jtpot{i}.variables).name); 
    sz(i)=length(jtpot{i}.variables); 
end

for s=1:S
    jtpotnum = whichpot(jtpot,sv(s),1); % find a single JT potential that contains the sympton
    margpot=sumpot(jtpot{jtpotnum},sv(s),0);
    p_symptom_jtree(s)=margpot.table(1);
end

% Print the results
fprintf(1, 'Junction tree algorithm\n');
for s=1:S
    fprintf(1,'p(s_%d=1) = %g\n',s,p_symptom_jtree(s)); 
end
disp(['Junction Tree ', num2str(t), ' seconds'])

%% More efficent algorithm
tstart=tic;
A = dag(pot);

for s=1:S
    margpot2 = sumpot(multpots([pot(sv(s)) pot(parents(A,sv(s)))]),sv(s),0);
    p_symptom_effalg(s)=margpot2.table(1);
end

t=toc(tstart);


fprintf(1, '\nEfficient algorithm\n');
for s=1:S
    fprintf(1,'p(s_%d=1) = %g\n',s, p_symptom_effalg(s)); 
end

[jtpot jtsep]=jtassignpot(setpot(pot,[sv(1:5) sv(6:10)],[1 1 1 1 1 2 2 2 2 2]),infostruct);
jtpot=absorption(jtpot,jtsep,infostruct); % do full round of absorption
fprintf(1, '\n');
disp('p(d_i = 1|s_1:10)')

for d=1:D
    jtpotnum = whichpot(jtpot,d,1) % find a single JT potential that contains the disease
    margpot=condpot(sumpot(jtpot(jtpotnum),d,0));
    p_disease_jtree(d) = margpot.table(1);
    fprintf(1,'p(d_%d = 1|s_1:10)=%g\n', d, p_disease_jtree(d));
end

disp(['Efficent inference takes ',num2str(t),' seconds.\n'])