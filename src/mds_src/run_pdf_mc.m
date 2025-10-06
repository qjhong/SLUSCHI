function [med] = run_pdf(system,avg_iter_vec,n_elms)

format shortG

for ii = 1:n_elms
for jj = 1:n_elms
    for j = 1:size(avg_iter_vec,2);
        avg_iter = round(avg_iter_vec(j))
        [S_xlogx,S_plogp,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s,elms] = pdf_mc(system,ii,jj,0,avg_iter,j);
        %[nstep/100,avg_iter,S_plogp]
        S_vec(j) = S_plogp;
    end
    close
    plot(avg_iter_vec,S_vec);hold on
    %for i = 1:size(S_vec,2)
    %    text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
    %end
    %n_int = size(S_vec,2)/4;
    %max_S = max(S_vec)+1.0;
    %med = zeros(4,1);
    %for j = 1:4
    %    i_start = 1+n_int*(j-1); i_end = n_int*(j);
    %    med(j) = median(S_vec(i_start:i_end));
    %    text(avg_iter_vec(i_start),max_S-0.2*j,strcat('median: ',num2str(med(j))));
    %end
    %min_S_vec = min(S_vec);
    ii,jj,S_vec
    ylim([0,S_vec(1)+0.2])
    title(strcat(replace(system,'_',' '),'K',' ',' S\_plogp',' ',replace(elms{ii},'_','\_'),' ',replace(elms{jj},'_','\_')))
    xlabel('n\_avg')
    ylabel('S [J/mol/K]')
    FS=12;
    set(gca,'FontSize',FS,'FontName','Times New Roman')%,'ytick',[1000,2000,3000,4000])
    set(findobj(gcf,'Type','text'),'FontSize',FS,'FontName','Times New Roman');
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Color', 'w');
    print('-f1','-r600','-dpng',strcat(system,'_S_plogp1_',elms{ii},'_',elms{jj}));
end
end
% if n_elms == 1
%     for j = 1:size(avg_iter_vec,2);
%         avg_iter = round(avg_iter_vec(j));
%         [S_xlogx,S_plogp,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s] = pdf_v5(system,nstep,1,1,0,avg_iter);
%         [nstep/100,avg_iter,S_plogp]
%         S_vec(j) = S_plogp;
%     end
%     close
%     plot(avg_iter_vec,S_vec);hold on
%     for i = 1:size(S_vec,2)
%         text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
%     end
%     n_int = size(S_vec,2)/4;
%     max_S = max(S_vec)+0.5;
%     for j = 1:4
%         i_start = 1+n_int*(j-1); i_end = n_int*(j);
%         med = median(S_vec(i_start:i_end));
%         text(avg_iter_vec(i_start+1),max_S-0.1*j,strcat('median: ',num2str(med)));
%     end
%     ylim([0,max_S+0.1])
%     print('-f1','-r600','-dpng',strcat(system,'_S_plogp1'));
% else if n_elms == 2
%     for j = 1:size(avg_iter_vec,2);
%         avg_iter = round(avg_iter_vec(j));
%         [S_xlogx,S_plogp12,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s] = pdf_v5(system,nstep,1,2,0,avg_iter);
%         [nstep/100,avg_iter,S_plogp12]
%         S_vec(j) = S_plogp12;
%     end
%     close
%     plot(avg_iter_vec,S_vec);hold on
%     for i = 1:size(S_vec,2)
%         text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
%     end
%     n_int = size(S_vec,2)/4;
%     max_S = max(S_vec)+0.5;
%     for j = 1:4
%         i_start = 1+n_int*(j-1); i_end = n_int*(j);
%         med = median(S_vec(i_start:i_end));
%         text(avg_iter_vec(i_start+1),max_S-0.1*j,strcat('median: ',num2str(med)));
%     end
%     ylim([0,max_S+0.1])
%     print('-f1','-r600','-dpng',strcat(system,'_S_plogp1_1_2'));
%     for j = 1:size(avg_iter_vec,2);
%         avg_iter = round(avg_iter_vec(j));
%         [S_xlogx,S_plogp21,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s] = pdf_v5(system,nstep,2,1,0,avg_iter);
%         [nstep/100,avg_iter,S_plogp21]
%         S_vec(j) = S_plogp21;
%     end
%     close
%     plot(avg_iter_vec,S_vec);hold on
%     for i = 1:size(S_vec,2)
%         text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
%     end
%     n_int = size(S_vec,2)/4;
%     max_S = max(S_vec)+0.5;
%     for j = 1:4
%         i_start = 1+n_int*(j-1); i_end = n_int*(j);
%         med = median(S_vec(i_start:i_end));
%         text(avg_iter_vec(i_start+1),max_S-0.1*j,strcat('median: ',num2str(med)));
%     end
%     ylim([0,max_S+0.1])
%     print('-f1','-r600','-dpng',strcat(system,'_S_plogp1_2_1'));
%     for j = 1:size(avg_iter_vec,2);
%         avg_iter = round(avg_iter_vec(j));
%         [S_xlogx,S_plogp11,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s] = pdf_v5(system,nstep,1,1,0,avg_iter);
%         [nstep/100,avg_iter,S_plogp11]
%         S_vec(j) = S_plogp11;
%     end
%     close
%     plot(avg_iter_vec,S_vec);hold on
%     for i = 1:size(S_vec,2)
%         text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
%     end
%     n_int = size(S_vec,2)/4;
%     max_S = max(S_vec)+0.5;
%     for j = 1:4
%         i_start = 1+n_int*(j-1); i_end = n_int*(j);
%         med = median(S_vec(i_start:i_end));
%         text(avg_iter_vec(i_start+1),max_S-0.1*j,strcat('median: ',num2str(med)));
%     end
%     ylim([0,max_S+0.1])
%     print('-f1','-r600','-dpng',strcat(system,'_S_plogp1_1_1'));
%     for j = 1:size(avg_iter_vec,2);
%         avg_iter = round(avg_iter_vec(j));
%         [S_xlogx,S_plogp22,S_plogp_merge1,S_plogp_merge2,R_cut,n_NN,avg_iter,R_x_s,R_anal_s,count_anal_s] = pdf_v5(system,nstep,2,2,0,avg_iter);
%         [nstep/100,avg_iter,S_plogp22]
%         S_vec(j) = S_plogp22;
%     end
%     close
%     plot(avg_iter_vec,S_vec);hold on
%     for i = 1:size(S_vec,2)
%         text(avg_iter_vec(i),S_vec(i),num2str(S_vec(i)))
%     end
%     n_int = size(S_vec,2)/4;
%     max_S = max(S_vec)+0.5;
%     for j = 1:4
%         i_start = 1+n_int*(j-1); i_end = n_int*(j);
%         med = median(S_vec(i_start:i_end));
%         text(avg_iter_vec(i_start+1),max_S-0.1*j,strcat('median: ',num2str(med)));
%     end
%     ylim([0,max_S+0.1])
%     print('-f1','-r600','-dpng',strcat(system,'_S_plogp1_2_2'));
% end

end
