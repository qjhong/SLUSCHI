addpath('replace_folder_here')
set(0, 'DefaultFigureVisible', 'off');
%nstep_vec = [3200];
%for i = 1:size(nstep_vec,2)
%nstep = nstep_vec(i);
system = ['replace_here'];
Nnu = 1000;
% step_unit = 2.4;
S_p_1_max = 0;
S_p_2_max = 0;
fid = fopen('vib.out','w');
[SS,S_p_1,S_p_2,F_1,F_2,nunu,S_id,n_int,natom,step_unit,n_elms,niter] = onephase_v6(system,Nnu);
S_p_1_max = max(S_p_1_max,S_p_1);
S_p_2_max = max(S_p_2_max,S_p_2);
fprintf(fid,'%f\n',S_p_1);
fprintf(fid,'%f\n',S_p_2);
fprintf(fid,'%f\n',S_id);
fprintf(fid,'%f\n',n_int);
fprintf(fid,'%f\n',natom);
S_p_1_final = S_p_1_max;
S_p_2_final = S_p_2_max;
S_p_1_final
S_p_2_final
S_id
n_int
natom
startExp = log10(round(10/step_unit)); % Starting value of a, in log10 scale
endExp = log10(round(min(2400/step_unit,niter/4))); % Ending value of b, in log10 scale
n = 20; % Number of points
avg_iter_vec = logspace(startExp, endExp, n);
%run_pdf_v6(system,avg_iter_vec,n_elms);
