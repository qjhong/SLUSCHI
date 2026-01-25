function [n_elms,n_atoms,mass,elms,n_atoms_total,stepsize,T,POS_LONG,LATT_LONG,VEL_LONG,step_unit] = read_files(filename)
    %%% READ IN FILES %%%
    fid = fopen(strcat('pos_',filename),'r');
    a = fscanf(fid,'%f');
    fclose(fid);
    n_values=size(a,1);
    
    fid = fopen(strcat('param_',filename),'r');
    data = fscanf(fid,'%f');
    n_elms = data(1); j = 1;
    n_atoms = zeros(n_elms,1);
    mass = zeros(n_elms,1);
    elms = {};
    for i = 1:n_elms
        j = j+1;
        n_atoms(i) = data(j);
    end
    for i = 1:n_elms
        j = j+1;
        mass(i) = data(j);
    end
    j = j+1; dt = data(j);
    j = j+1; n_atoms_total = data(j);
    for i = 1:n_elms
        elms{i} = fgetl(fid);
    end
    fclose(fid);
    
    niter=n_values/6/n_atoms_total;
    fid = fopen(strcat('latt_',filename),'r');
    b = fscanf(fid,'%f');
    fclose(fid);
    
    fid = fopen(strcat('step_',filename),'r');
    stepsize = fscanf(fid,'%f');
    unique_stepsize = unique(stepsize);
    A = round(unique_stepsize*10);
    if size(A,1) == 1
        gcd_step = A;
    else
        gcd_step = gcd(A(1),A(2));
        for i = 3:size(A,1)
            gcd_step = gcd(gcd_step,A(i));
        end
    end
    step_unit = gcd_step/10;
    fclose(fid);
    
    mass = mass/1000/6.023e23;
    % dt = dt*1e-15; % 3fs
    output = strsplit(filename, "_");
    T=str2num(output{2});
    
    %niter
    POS=zeros(3,n_atoms_total,niter); % in A
    kk=0;
    for iter=1:niter
        for iatom=1:n_atoms_total
            for j=1:6
                kk=kk+1;
                if (j<=3)
                POS(j,iatom,iter) = a(kk);
                end
            end
        end
    end
    
    LATT=zeros(3,3,niter/80); % in A
    kk=0;
    for iter=1:niter/80
        for iatom=1:3
            for j=1:6
                kk=kk+1;
                if (j<=3)
                LATT(j,iatom,iter) = b(kk);
                end
            end
        end
    end
    %size(stepsize)

    % Calculate Velocity
    VEL=zeros(3,n_atoms_total,niter); % in A/3 fs
    for iter=1:niter-1
        VEL(:,:,iter) = POS(:,:,iter+1)-POS(:,:,iter);
        for iatom=1:n_atoms_total
            step_cur = stepsize( floor((iter)/80) + 1 );
            if norm(VEL(:,iatom,iter)) > 1.
                stopflag = 0;min_norm = 100;
                for i=-1:1
                for j=-1:1
                for k=-1:1
                    tmp = VEL(:,iatom,iter);
                    tmp = tmp + LATT( :, :, floor((iter-1)/80) + 1 ) * [i;j;k];
                    tmp_norm = norm(tmp);
                    min_norm = min(min_norm,tmp_norm);
                    if norm(tmp) < 1.
                        stopflag = 1;
                        break
                    end
                end
                if stopflag==1
                    break
                end
                end
                if stopflag==1 
                    break
                end
                end
                if stopflag == 0
                    VEL(:,iatom,iter)
                    LATT( :, :, floor((iter-1)/80) + 1 )
                    min_norm
                end
                VEL(:,iatom,iter) = tmp;
            end
        end
    end
    VEL = VEL *1e-10/(step_cur*1e-15);
    
    % VEL_v2=zeros(3*n_atoms_total,niter); % in A/3 fs
    % for i=1:3
    % for j=1:n_atoms_total
    % for k=1:niter
    %     VEL_v2(i+3*(j-1),k) = VEL(i,j,k);
    % end
    % end
    % end
    for k=80:80:niter-1
        VEL(:,:,k) = (VEL(:,:,k-1) + VEL(:,:,k+1) )/2;
    end

    %%% CHECK %%%
    if size(a,1)/6/n_atoms_total ~= size(b,1)/18*80
        error('Dimension Error: POS and LATT');
    end
    if size(b,1)/18 ~= size(stepsize,1)
        error('Dimension Error: LATT and STEP');
    end

    %%% NOW SET UP POS AND LATT, ASSUMING STEPSIZE NOT UNIFORM %%%
    niter_long = round(sum(stepsize)/step_unit*80);
    POS_LONG = zeros(3,n_atoms_total,1);
    VEL_LONG = zeros(3,n_atoms_total,1);
    LATT_LONG = zeros(3,3,1);
    POS_LONG(:,:,1) = POS(:,:,1);
    VEL_LONG(:,:,1) = VEL(:,:,1);
    LATT_LONG(:,:,1) = LATT(:,:,1);
    iter_pos_long = 1;
    for i = 2:size(POS,3)
        step_cur = stepsize( floor((i-1)/80) + 1 );
        dPOS = POS(:,:,i) - POS(:,:,i-1); 
        dVEL = VEL(:,:,i) - VEL(:,:,i-1); 
        latt_cur = LATT(:,:,floor((i-1)/80) + 1);
        dPOS = convert_in_cell(dPOS,latt_cur);
        for j = 1:round(step_cur/step_unit)
            iter_pos_long = iter_pos_long + 1;
            POS_LONG(:,:,iter_pos_long) = POS(:,:,i) + dPOS / round(step_cur/step_unit) * j;
            VEL_LONG(:,:,iter_pos_long) = VEL(:,:,i) + dVEL / round(step_cur/step_unit) * j;
            LATT_LONG(:,:,iter_pos_long) = latt_cur;
        end
    end
end

function VEL_new = convert_in_cell(VEL,LATT)
    VEL_new = VEL;
    for iatom=1:size(VEL,2)
    if norm(VEL(:,iatom)) > 1.
        stopflag = 0;min_norm = 100;
        for i=-1:1
        for j=-1:1
        for k=-1:1
            tmp = VEL(:,iatom);
            tmp = tmp + LATT * [i;j;k];
            tmp_norm = norm(tmp);
            if tmp_norm < min_norm
                min_norm = min(min_norm,tmp_norm);
                min_tmp = tmp;
            end
            if norm(tmp) < 5.
                stopflag = 1;
                break
            end
        end
        if stopflag==1
            break
        end
        end
        if stopflag==1 
            break
        end
        end
        if stopflag == 0
            %VEL(:,iatom)
            %LATT( :, :, floor((iter-1)/80) + 1 )
            %min_norm
        end
        VEL_new(:,iatom) = min_tmp;
    end
    end
end
