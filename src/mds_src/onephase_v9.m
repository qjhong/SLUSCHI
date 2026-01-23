function [SS,S_p_1,S_p_2,F_1,F_2,nunu,S_id,n_int,n_atoms_total,step_unit,n_elms,niter] = onephase_v6(filename,Nnu)

h=6.626e-34;
eV = 1.602e-19;
k_B=1.38e-23;

[n_elms,n_atoms,mass,elms,n_atoms_total,stepsize,T,POS,LATT,VEL,step_unit] = read_files(filename);
step_unit
niter = size(POS,3);
dt=step_unit*1e-15;

VEL_v2=zeros(3*n_atoms_total,niter); % in A/3 fs
for i=1:3
for j=1:n_atoms_total
for k=1:niter
    VEL_v2(i+3*(j-1),k) = VEL(i,j,k);
end
end
end

nsteps = floor(niter/2)-100;
nsteps = floor(nsteps/2)*2;
% nsteps = 10000;
t0 = nsteps; tau = nsteps/2;
% Nnu=200;
nu_max=10e13;
NU = linspace(0,nu_max,Nnu);
S = zeros(Nnu,n_elms);

idx_start = 1;
for i_elm = 1:n_elms
    idx_end = idx_start - 1 + n_atoms(i_elm)*3;
    VEL_v3 = VEL_v2(idx_start:idx_end,:);
    % VEL_MAT = zeros(2*t0,2*t0);
    % for Iiter = t0-tau:t0+tau-1
    % for Jiter = Iiter-tau+1:Iiter+tau
    %     VEL_MAT(Iiter,Jiter) = VEL_v3(:,Iiter)'*VEL_v3(:,Jiter);
    % end
    % VEL_MAT(Iiter,Iiter-tau+1:Iiter+tau) = VEL_v3(:,Iiter)'*VEL_v3(:,Iiter-tau+1:Iiter+tau);
    % end    
    % for inu = 1:Nnu
        % nu = NU(inu);
        factor = -2*pi*NU'*dt;
        for Iiter = t0-tau:t0+tau-1
        % for Jiter = Iiter-tau+1:Iiter+tau
        %     S(inu,i_elm) = S(inu,i_elm) + VEL_MAT(Iiter,Jiter)*dt^2*cos(-2*pi*nu*(Iiter-Jiter)*dt);
        % end
            Jiter = Iiter-tau+1:Iiter+tau;
            COSJ = cos(factor*(Iiter-Jiter));
            % S(inu,i_elm) = S(inu,i_elm) + VEL_MAT(Iiter,Iiter-tau+1:Iiter+tau)*COSJ';
            S(:,i_elm) = S(:,i_elm) + COSJ*(VEL_v3(:,Iiter)'*VEL_v3(:,Iiter-tau+1:Iiter+tau))';
        end
        S(:,i_elm) = S(:,i_elm)*dt^2;
    % end
    S(:,i_elm) = S(:,i_elm)/2/(tau*dt)*mass(i_elm)*2/k_B/T;
    idx_start = idx_start + n_atoms(i_elm)*3;
end
% VEL_v2_fft = fft(VEL_v2');
% sumS = zeros(size(VEL_v2_fft,1),1);
% for i = 1:natom*3
%     sumS = sumS + VEL_v2_fft(:,i).*conj(VEL_v2_fft(:,i));
% end

V=0;
for i = 1:size(LATT,3)
    V = V + det(LATT( :, :, i));
end
V = V/(size(LATT,3));

sumS = S;

S_p_1 = zeros(n_elms,1);
S_p_2 = zeros(n_elms,1);
S_id = zeros(n_elms,1);
n_int = zeros(n_elms,1);
for i_elm = 1:n_elms
    plot(sumS(:,i_elm));hold on
    n_int(i_elm) = sum(sumS(:,i_elm))*nu_max/(Nnu-1);
end

n = n_atoms_total;
V = V/n*1e-30;
fid = fopen('phonon_dos_v9.out','w');
for i_elm = 1:n_elms
lambda = sqrt(h^2/(2*pi*mass(i_elm)*k_B*T));
S_id(i_elm)=(5/2 + log(V/lambda^3) )*k_B*6.023e23;
%phonon
nn = 10000;
for i = 1:1
%nn = nn*10;
nunu = linspace(0,nu_max,nn);
%sumS(:,i_elm)
SS = spline(NU,sumS(:,i_elm),nunu);
F_1 = zeros(nn,1);
F_2 = zeros(nn,1);
E_1 = zeros(nn,1);
area = 0;
max_SS=0;idx_SS=0;
flag_correction=1;
for i=1:nn-1
    %nu_curr = ( nu(i)+nu(i+1) ) / 2;
    %F(i+1) = F(i) + k * T * log( 2 * sinh( h*nu_curr / (2*k*T) ) ) * (S(i)+S(i+1))/2*nu_max/(n-1);
    %F(i+1) = k * T * log( 2 * sinh( h*nu_curr / (2*k*T) ) ) * (S(i)+S(i+1))/2*nu_max/(n-1);
    nu_curr = ( nunu(i)+nunu(i+1) ) / 2;
    F_1(i+1) = T * k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) )  * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    F_2(i+1) = T * max(-S_id(i_elm)/3/6.023e23+k_B, k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) ) ) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    E_1(i+1) = h*nu_curr/2 *  (exp(h*nu_curr/k_B/T)+1) / (exp(h*nu_curr/k_B/T)-1) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    area = area + (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    % 3 translational degrees of freedom are missing in VASP
    % assuming: 1. they are the softest modes in DOS
    %           2. they are equally distributed among elements
    % 3 DOF distributed in n_atoms(i_elm) / sum(n_atoms)
    % area final value is 3 * n_atoms(i_elm)
    if flag_correction==1 && area > 3 * n_atoms(i_elm) / sum(n_atoms) % time to keep a copy of correction
	area_c = area;
	E_1_c = sum(E_1(1:i+1));
	F_1_c = sum(F_1(1:i+1));
	F_2_c = sum(F_2(1:i+1));
	flag_correction=0;
    end
    if SS(i) > max_SS
        max_SS = SS(i);
        idx_SS = i;
    end
end
area_id = 0;
for i=1:idx_SS
    area_id = area_id + (SS(i)+SS(i+1))/2*nu_max/(nn-1);
end
close
plot(nunu,SS)
grid on
FS=16;
xlabel('$\nu$ [s$^{-1}$]','Interpreter','Latex')
ylabel('$DoS$','Interpreter','Latex')
title(strcat(replace(filename,'_','\_'),'K ','\_',replace(elms(i_elm),'_','\_')))
set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r600','-dpng',strcat(filename,'_phonon_DOS_',elms{i_elm}));
close
plot(nunu,k_B*(1 - log( 2 * sinh( h*nunu / (2*k_B*T) ) ))*6.023e23)
%pause(10)
close
plot(nunu,F_1/area*3/eV,'b')
hold on
plot(nunu,F_2/area*3/eV,'r')
grid on
xlabel('$\nu$ [s$^{-1}$]','Interpreter','Latex')
ylabel('$dF_{vib}/d\nu$','Interpreter','Latex')
title(strcat(replace(filename,'_','\_'),'K ','\_',replace(elms(i_elm),'_','\_')))
set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r600','-dpng',strcat(filename,'_S_',elms{i_elm}));

area0 = area;
area = 0;
int0 = zeros(nn,1);
int1 = zeros(nn,1);
for i=1:nn-1
    nu_curr = ( nunu(i)+nunu(i+1) ) / 2;
    area = area + (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    % dos = k * x^2;
    % int = 1/3 * k * x^3;
    k = area*3/(nunu(i+1)^3);
    for j = 1:i
        nu_curr = ( nunu(j)+nunu(j+1) ) / 2;
        SS_curr = k * nu_curr^2;
        F_1_corr(j+1) = T * k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) )  * SS_curr *nu_max/(nn-1);
        E_1_corr(j+1) = h*nu_curr/2 *  (exp(h*nu_curr/k_B/T)+1) / (exp(h*nu_curr/k_B/T)-1) * SS_curr *nu_max/(nn-1);
    end
    for j = i+1:nn-1
        nu_curr = ( nunu(j)+nunu(j+1) ) / 2;
        F_1_corr(j+1) = T * k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) )  * (SS(j)+SS(j+1))/2*nu_max/(nn-1);
        E_1_corr(j+1) = h*nu_curr/2 *  (exp(h*nu_curr/k_B/T)+1) / (exp(h*nu_curr/k_B/T)-1) * (SS(j)+SS(j+1))/2*nu_max/(nn-1);
    end
    int0(i) = -(sum(F_1)/area0*3/eV-sum(E_1)/area0*3/eV)/T*96485;
    int1(i) = -(sum(F_1_corr)/area0*3/eV-sum(E_1_corr)/area0*3/eV)/T*96485;
    % F_1(i+1) = T * k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) )  * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    % F_2(i+1) = T * max(-S_id(i_elm)/3/6.023e23+k_B, k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) ) ) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    % E_1(i+1) = h*nu_curr/2 *  (exp(h*nu_curr/k_B/T)+1) / (exp(h*nu_curr/k_B/T)-1) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
end
close
plot(SS)
pause(1)
close
plot(int0)
hold on
plot(int1)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r600','-dpng',strcat(filename,'_S_corr_',elms{i_elm}));
[v,idx]=max(int1)
area = 0;
for i=1:idx
    nu_curr = ( nunu(i)+nunu(i+1) ) / 2;
    area = area + (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    k = area*3/(nunu(i+1)^3);
end
for j = 1:idx
    nu_curr =  nunu(j);
    SS(j) = k * nu_curr^2;
end
close
plot(nunu,SS)
grid on
FS=16;
xlabel('$\nu$ [s$^{-1}$]','Interpreter','Latex')
ylabel('$DoS$','Interpreter','Latex')
title(strcat(replace(filename,'_','\_'),'K ','\_',replace(elms(i_elm),'_','\_')))
set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
print('-f1','-r600','-dpng',strcat(filename,'_phonon_DOS_corr_',elms{i_elm}));

area = 0;
flag_correction=1;
for i=1:nn-1
    %nu_curr = ( nu(i)+nu(i+1) ) / 2;
    %F(i+1) = F(i) + k * T * log( 2 * sinh( h*nu_curr / (2*k*T) ) ) * (S(i)+S(i+1))/2*nu_max/(n-1);
    %F(i+1) = k * T * log( 2 * sinh( h*nu_curr / (2*k*T) ) ) * (S(i)+S(i+1))/2*nu_max/(n-1);
    nu_curr = ( nunu(i)+nunu(i+1) ) / 2;
    F_1(i+1) = T * k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) )  * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    F_2(i+1) = T * max(-S_id(i_elm)/3/6.023e23+k_B, k_B * log( 2 * sinh( h*nu_curr / (2*k_B*T) ) ) ) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    E_1(i+1) = h*nu_curr/2 *  (exp(h*nu_curr/k_B/T)+1) / (exp(h*nu_curr/k_B/T)-1) * (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    area = area + (SS(i)+SS(i+1))/2*nu_max/(nn-1);
    % 3 translational degrees of freedom are missing in VASP
    % assuming: 1. they are the softest modes in DOS
    %           2. they are equally distributed among elements
    % 3 DOF distributed in n_atoms(i_elm) / sum(n_atoms)
    % area final value is 3 * n_atoms(i_elm)
    if flag_correction==1 && area > 3 * n_atoms(i_elm) / sum(n_atoms) % time to keep a copy of correction
	area_c = area;
	E_1_c = sum(E_1(1:i+1));
	F_1_c = sum(F_1(1:i+1));
	F_2_c = sum(F_2(1:i+1));
	flag_correction=0;
    end
    if SS(i) > max_SS
        max_SS = SS(i);
        idx_SS = i;
    end
end

E = 3*k_B*T/eV
%E_1 = (sum(E_1))/area*3/eV
E_1 = (sum(E_1)+E_1_c)/(area+area_c)*3/eV
%FF_1 = (sum(F_1)-S_id/6.023e23*T)/area*3/eV;
%FF_1 = (sum(F_1))/area*3/eV;
%S_p_1(i_elm) = -(FF_1-E_1)/T*96485
FF_1 = (sum(F_1)+F_1_c)/(area+area_c)*3/eV;
%S_p_1(i_elm) = -(FF_1-E)/T*96485;
S_p_1(i_elm) = -(FF_1-E_1)/T*96485;
%FF_2 = (sum(F_2)-S_id/6.023e23*T)/area*3/eV;
%FF_2 = (sum(F_2))/(area)*3/eV;
%S_p_2(i_elm) = -(FF_2-E_1)/T*96485
FF_2 = (sum(F_2)+F_1_c)/(area+area_c)*3/eV;
%S_p_2(i_elm) = -(FF_2-E)/T*96485;
S_p_2(i_elm) = -(FF_2-E_1)/T*96485;
%area_id / area
fprintf(fid, '%6.2e\n', nunu);
fprintf(fid, '%6.2e\n', SS);
fprintf(fid, '%6.2e\n', F_1/area*3/eV);
fprintf(fid, '%6.2e\n', F_2/area*3/eV);
end
end
fclose(fid);
