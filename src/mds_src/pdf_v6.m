function [S_xlogx,S_plogp,S_plogp_merge1,S_plogp_merge2,R_cut0,n_NN,avg_iter,R_x,R_anal,count_anal,elms] = pdf_v6(filename,idx_in1,idx_in2,flag_solid,avg_iter,j_iter)

h=6.626e-34;
eV = 1.602e-19;
k_B=1.38e-23;

[n_elms,n_atoms,mass,elms,n_atoms_total,stepsize,T,POS,LATT,VEL,step_unit] = read_files(filename);
niter = size(POS,3);

% === Downsampling controls (tune these) ===
% Process every FRAME_STRIDE-th averaged frame (1 = all, 2 = half, 5 = 1/5, etc.)
FRAME_STRIDE = max(1,round(niter/avg_iter/300))

% Calculate velocity and position
R_max = 10;
n_R = 100;
dR = R_max/n_R;
R_x = linspace(dR,R_max,n_R);
R_anal = zeros(n_R,n_elms,n_elms);
%avg_iter = 160;
n_atoms_acc = n_atoms;
for i = 2:size(n_atoms,1)
    n_atoms_acc(i) = n_atoms_acc(i-1) + n_atoms_acc(i);
end
for iter=1:niter
    if mod(iter,avg_iter)==1
        POS0 = POS(:,:,iter);
    else
        DIST(:,:) = POS(:,:,iter)-POS0;
        latt = LATT( :, :, iter );
        % minimal-image mapping using fractional coordinates + round()
        % R_frac = inv(latt) * DIST  gives fractional coords (3 x n_atoms)
        invlatt = inv(latt);                   % (3x3)
        R_frac_all = invlatt * DIST;           % (3 x n_atoms)
	% nearest integer lattice vector for each atom
   	shift_cells = round(R_frac_all);        % (3 x n_atoms) of integers

   	% convert integer shift back to Cartesian displacement and subtract
   	% DIST_new = DIST - latt * shift_cells
   	DIST = DIST - latt * shift_cells;      % (3 x n_atoms)

   	% update POS for this iter
   	POS(:,:,iter) = POS(:,:,iter-1) + POS0 + DIST;
        %for iatom=1:n_atoms_total
        %if norm(DIST(:,iatom)) > 1.
        %    stopflag = 0;min_norm = 100;
        %    for i=-1:1
        %    for j=-1:1
        %    for k=-1:1
        %        tmp = DIST(:,iatom);
        %        % tmp = tmp + LATT( :, :, floor((iter-1)/80) + 1 ) * [i;j;k];
        %        tmp = tmp + latt * [i;j;k];
        %        tmp_norm = norm(tmp);
        %        if tmp_norm < min_norm
        %            min_norm = min(min_norm,tmp_norm);
        %            min_tmp = tmp;
        %        end
        %        if norm(tmp) < 5.
        %            stopflag = 1;
        %            break
        %        end
        %    end
        %    if stopflag==1
        %        break
        %    end
        %    end
        %    if stopflag==1 
        %        break
        %    end
        %    end
        %    if stopflag == 0
        %        %VEL(:,iatom)
        %        %LATT( :, :, floor((iter-1)/80) + 1 )
        %        %min_norm
        %    end
        %    DIST(:,iatom) = min_tmp;
        %end
        %end
        %POS(:,:,iter) = POS(:,:,iter-1) + POS0 + DIST;
    end
    %calculate radial distribution
    if mod(iter,avg_iter*FRAME_STRIDE) == 0
        POS(:,:,iter) = POS(:,:,iter)/avg_iter;
        latt = LATT( :, :, iter );
        invlatt = inv(latt);
        POS_iter = POS(:,:,iter);
        for iatom=1:n_atoms_total
        for jatom=iatom+1:n_atoms_total
            R = POS_iter(:,iatom)-POS_iter(:,jatom);
            R_frac = invlatt * R;
            %for i = 1:3
            %    while R_frac(i) > 0.5
            %        R_frac(i) = R_frac(i) - 1.;
            %    end
            %    while R_frac(i) < -0.5
            %        R_frac(i) = R_frac(i) + 1.;
            %    end
            %end
	    %for i1=-1:1
	    %for i2=-1:1
	    %for i3=-1:1
	    % old: loop for i =1:3 adjust >0.5 or < -0.5
            R_frac = R_frac - round(R_frac);   % now components in [-0.5,0.5]
            R = latt * R_frac;
            min_norm = norm(R);
            if min_norm < R_max
                idx = ceil(min_norm/dR);
                for i_elm = 1:n_elms
                    if iatom <= n_atoms_acc(i_elm)
                        break
                    end
                end
                for j_elm = 1:n_elms
                    if jatom <= n_atoms_acc(j_elm)
                        break
                    end
                end                
	        R_anal(idx,i_elm,j_elm) = R_anal(idx,i_elm,j_elm) + 1;
            end
	    %end
	    %end
	    %end
        end
        end
    end
    %iter
end

V=0;
for i = 1:niter/80
    V = V + det(LATT( :, :, i));
end
V = V/(niter/80);

close
%plot(R_x,R_anal)
%pause(10)
for i_elm = 1:n_elms
for j_elm = i_elm:n_elms
for i = 1:n_R
    R_anal(i,i_elm,j_elm) = R_anal(i,i_elm,j_elm) / (4*pi*R_x(i)^2*dR * (n_atoms(i_elm)-1)/V) / (niter/avg_iter) / n_atoms(j_elm);
end
if i_elm == j_elm
    R_anal(:,i_elm,j_elm) = R_anal(:,i_elm,j_elm) * 2.;
end
end
end
close
res_rr=R_x';
for i_elm = 1:n_elms
for j_elm = i_elm:n_elms
    res_rr=[res_rr, R_anal(:,i_elm,j_elm)];
    plot(R_x,R_anal(:,i_elm,j_elm));hold on
end
end
FS=12;
xlabel('$r$ [\AA]','Interpreter','Latex')
ylabel('$g(r)$','Interpreter','Latex')
title(strcat(replace(filename,'_',' '),'K'))
set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r300','-dpng',strcat(filename,'_R_anal','_',num2str(j_iter)));

% define cutoff
flag_increase = 0;
flag_decrease = 0;
n_NN = 0;
R_anal0 = R_anal;
R_anal = R_anal0(:,min(idx_in1,idx_in2),max(idx_in1,idx_in2));
R_max = max(R_anal(1:n_R/2));
%R_anal
for i = 1:n_R
    %if R_anal(i)>R_anal(i-1)+1e-10 && flag_increase == 0
    if R_anal(i)>R_max-1e-10 && flag_increase == 0
        flag_increase = 1;
    end
    if i>1 && R_anal(i)<R_anal(i-1)-1e-10 && flag_increase == 1 && R_anal(i-1) > 0.1 * R_max % make sure this is a large peak
        flag_decrease = 1;
        %i
        %[R_anal(i-1:i)]
    end
    %if flag_increase == 1 && flag_decrease == 1 && R_anal(i)>R_anal(i-1)+1e-10 && max(R_anal(i-1:i)-R_anal(i:i+1)) < 0. % for two "merging" peaks which approach zero after R_cut
    if i>1 && flag_increase == 1 && flag_decrease == 1 && R_anal(i)>R_anal(i-1)+1e-10 && max(R_anal(i-1:i)-R_anal(i:i+1)) < 0. % for two "merging" peaks which approach zero after R_cut
      %R_anal(i-1:i)-R_anal(i:i+1)
        R_cut = R_x(i-1);
        if R_cut > 5
	  R_cut
          R_anal
        end
        for j = 1:i-1
            n_NN = n_NN + 4*pi*R_x(j)^2*dR * (n_atoms_total-1)/V * R_anal(j);
        end
        break
    end
end
n_NN;
R_cut0 = R_cut;

% calcualte radial distribution
res_count=[ ];
n_R_cut = 0;
S_R_cut = zeros(n_R_cut,1);
for i_R_cut = 0:n_R_cut
    R_cut = R_cut - 0.05;
n_max = 100;
count_anal = zeros(n_max,1);
%n_atoms_acc(idx_in1)-n_atoms(idx_in1)+1,n_atoms_acc(idx_in1)
%n_atoms_acc(idx_in2)-n_atoms(idx_in2)+1,n_atoms_acc(idx_in2)
for iter=avg_iter*FRAME_STRIDE:avg_iter*FRAME_STRIDE:niter
    latt = LATT( :, :, iter );
    invlatt = inv(latt);
    POS_iter = POS(:,:,iter);
    for iatom=n_atoms_acc(idx_in1)-n_atoms(idx_in1)+1:n_atoms_acc(idx_in1)
        count = 0;
    R_list = ones(1,1000)*100;
    for jatom=n_atoms_acc(idx_in2)-n_atoms(idx_in2)+1:n_atoms_acc(idx_in2)
        if iatom ~= jatom
        R = POS_iter(:,iatom)-POS_iter(:,jatom);
	% minimal-image mapping (vectorized, puts fractional coords in [-0.5,0.5])
        R_frac = invlatt * R;
        %R_frac = R_frac - round(R_frac);
        %R_frac = invlatt * R;
        %for i = 1:3
        %    while R_frac(i) > 0.5
        %        R_frac(i) = R_frac(i) - 1.;
        %    end
        %    while R_frac(i) < -0.5
        %        R_frac(i) = R_frac(i) + 1.;
        %    end
        %end
	%for i1=-1:1
	%for i2=-1:1
	%for i3=-1:1
	% old: loop for i =1:3 adjust >0.5 or < -0.5
        R_frac = R_frac - round(R_frac);   % now components in [-0.5,0.5]
        R = latt * R_frac ;
        min_norm = norm(R);
        %R_list(jatom) = min_norm;
        if min_norm < R_cut
            count = count + 1;
            R_list(count) = min_norm;
        end
        end
        %end
        %end
        %end
    end
    % % for solid, adjust count
    % if flag_solid == 1
    %     % [count,R_cut]
    %     R_list = sort(R_list);
    %     %R_list(1:10)
    %     %    iter,R_cut,count,R_list(1:14)
    %     d0 = R_list(count+1) - R_list(count);
    %     d1 = R_list(count) - R_list(count-1);
    %     d2 = R_list(count+2) - R_list(count+1);
    %     if d1 > d0
    %         count = count - 1;
    %     else 
    %         if d2 > d0
    %             count = count + 1;
    %         end
    %     end
    % end
    if idx_in1==3 & idx_in2 ==3 & count ~= 12
	%R_list(1:14)
	%quit
        %R_list = sort(R_list);
        %iter,R_cut,count,R_list(1:14)
    end
    %count
    if count > 0
    count_anal(count) = count_anal(count) + 1;
    %[iatom,count]
    end
    end
    %iter
end
close

S = 0;
count_anal;
count_anal = count_anal / sum(count_anal);
count_anal_max = max(count_anal);
%plot(count_anal)
%print('-f1','-r600','-dpng',strcat(filename,'_count_anal'));

for i = 1:n_max
    x = i / max(n_NN,i+1);
    S = S + count_anal(i) * (x*log(x) + (1-x)*log(1-x)) / x;
end
S;
S_xlogx = -S*8.314;
S_plogp = sum(count_anal.*log(count_anal+0.00001))*8.314/2*-1;
S_R_cut(i_R_cut+1) = S_plogp;
end
S_R_cut;
[S_plogp,idx] = min(S_R_cut);
R_cut = R_cut0 - 0.05*(idx);
n_max = 100;
% count_anal = zeros(n_max,1);
% for iter=avg_iter:avg_iter:niter
%     for iatom=n_atoms_acc(idx_in1)-n_atoms(idx_in1)+1:n_atoms_acc(idx_in1)
%         count = 0;
%     R_list = ones(1,1000)*100;
%     for jatom=n_atoms_acc(idx_in2)-n_atoms(idx_in2)+1:n_atoms_acc(idx_in2)
%         if iatom ~= jatom
%         R = POS(:,iatom,iter)-POS(:,jatom,iter);
%         latt = LATT( :, :, floor((iter-1)/80) + 1 );
%         R_frac = inv(latt) * R;
%         for i = 1:3
%             while R_frac(i) > 0.5
%                 R_frac(i) = R_frac(i) - 1.;
%             end
%             while R_frac(i) < -0.5
%                 R_frac(i) = R_frac(i) + 1.;
%             end
%         end
%         R = latt * R_frac;
%         min_norm = norm(R);
%         R_list(jatom) = min_norm;
%         if min_norm < R_cut
%             count = count + 1;
%         end
%         end
%     end
%     % % for solid, adjust count
%     % if flag_solid == 1
%     %     % [count,R_cut]
%     %     R_list = sort(R_list);
%     %     %R_list(1:10)
%     %     %    iter,R_cut,count,R_list(1:14)
%     %     d0 = R_list(count+1) - R_list(count);
%     %     d1 = R_list(count) - R_list(count-1);
%     %     d2 = R_list(count+2) - R_list(count+1);
%     %     if d1 > d0
%     %         count = count - 1;
%     %     else 
%     %         if d2 > d0
%     %             count = count + 1;
%     %         end
%     %     end
%     % end
%     % if count > 0
%     %     %count
%     % count_anal(count) = count_anal(count) + 1;
%     % end
%     end
%     %iter
% end
close

S = 0;
count_anal;
count_anal = count_anal / sum(count_anal);
bar(count_anal)
res_count=[res_count, count_anal];
xlabel('$n$','Interpreter','Latex')
ylabel('$P_n$','Interpreter','Latex')
ylim([0,1])
title(strcat(replace(filename,'_',' '),'K ','\_',replace(elms{idx_in1},'_','\_'),'\_',replace(elms{idx_in2},'_','\_')))
set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r300','-dpng',strcat(filename,'_count_anal_',elms{idx_in1},'_',elms{idx_in2},'_',num2str(j_iter)));

for i = 1:n_max
    x = i / max(n_NN,i+1);
    S = S + count_anal(i) * (x*log(x) + (1-x)*log(1-x)) / x;
end
S;
S_xlogx = -S*8.314;
S_plogp = sum(count_anal.*log(count_anal+0.00001))*8.314/2*-1;
% sigma=1;[tmp,nu]=max(count_anal);
% n = linspace(1,n_max,n_max);
% f = 1 ;
% vec = [f, sigma, nu]; vec0 = vec;
% diff = vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 ) - count_anal';
% cost0 = diff*diff';
% n_opt = 1000;
% diff_size = 0.01;
% step_size = 0.001;
% for i = 1:n_opt
%     % gradient
%     grad = zeros(1,3);
%     for j = 1:3
%         vec = vec0;
%         vec(j) = vec(j) + diff_size;
%         diff = vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 ) - count_anal';
%         cost = diff*diff';
%         grad(j) = (cost - cost0)/diff_size;
%     end
%     vec = vec0 - grad * step_size;
%     diff = vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 ) - count_anal';
%     cost0 = diff*diff';
%     vec0 = vec;
%     close
%     plot(count_anal)
%     hold on
%     plot(vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 ))
% end
% count_anal2 = vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 );
% S_plogp_fit = sum(count_anal2.*log(count_anal2+0.00001))*8.314/2*-1;

% vec(2) = vec(2) - 0.5;
% plot(vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 ))
% count_anal3 = vec(1) / vec(2) / sqrt(2*pi) * exp( -1/2 * ((n-vec(3))/vec(2)).^2 );
% S_plogp_fit_narrow = sum(count_anal3.*log(count_anal3+0.00001))*8.314/2*-1;

count_anal4 = zeros(size(count_anal));
for i = 1:size(count_anal,1)
    if mod(i,2)==0
        count_anal4(i) = count_anal(i) + count_anal(i-1);
    end
end
S_plogp_merge1 = sum(count_anal4.*log(count_anal4+0.00001))*8.314/2*-1;

count_anal5 = zeros(size(count_anal));
for i = 1:size(count_anal,1)
    if mod(i,2)==1 && i>1
        count_anal5(i) = count_anal(i) + count_anal(i-1);
    end
end
S_plogp_merge2 = sum(count_anal5.*log(count_anal5+0.00001))*8.314/2*-1;
%[S_xlogx,S_plogp,S_plogp_merge1,S_plogp_merge2]

isave=99;
if isave > 0;
   ssname = strcat(filename,'_count_', elms{idx_in1},'_', elms{idx_in2}, '.txt');
   eval(['save ' ssname ' res_count -ascii']);

   ssname=strcat(filename,'_RDF.txt'); eval(['save ' ssname ' res_rr -ascii']);
end
