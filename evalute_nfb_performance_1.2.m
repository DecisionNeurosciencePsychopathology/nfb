function evalute_nfb_performance(out)
%Function will evalute nfb subjects responses, place resutls into graphics
%for ease of viewing and assisting with evaluation of subejcts perofmrance 
%on a case by case basis as well as experiment wide performance.

%Grab field names -- we only want to use Master Data Files
fnames = fieldnames(out);

for i = 1:length(fnames)
   if strfind(fnames{i},'MDF')>0 
       %MDF data
       data=out.(fnames{i});
       
       %Expermient wide plot
       plot_behavior(data);
       
       %TODO
       %Grab the indivual subjects here (regrep(Participant numbers, unqiue
       %it?, cellfun(@(x) strcmp...) feed to
       %bar_plot_behavior(data,graph_options)
       
       %Individual subject loop
       subjs=regexp(out.(fnames{i}).Participant,'\d{3}','match');
       subjs=unique([subjs{:}]');
       for j = 1:length(subjs)
           subjs_data_idx = cellfun(@(x) ~isempty(strfind(x,subjs{j})),data.Participant,'UniformOutput',0);
           subjs_data = data(cell2mat(subjs_data_idx),:);
           plot_behavior(subjs_data);
       end
       
   end
end


function plot_behavior(data)
%TODO
%Make this the generic plot funciton with a second arguement graph_options
%graph_options will have things like Infusion_title, subject(s)_name,
%not_infusion_title, type_of_plot(bar line ... )ect...
infusion_trials = cellfun(@(x) strcmpi(x,'A') | strcmpi(x,'C') ,data.Infusion);
will_imp_resp_inf = data.WillImpRespNum(infusion_trials);
will_imp_resp_no_inf = data.WillImpRespNum(~infusion_trials);
figure(1)
subplot(1,2,1)
bar([sum(will_imp_resp_inf==1),sum(will_imp_resp_inf==0), sum(isnan(will_imp_resp_inf))])
subplot(1,2,2)
bar([sum(will_imp_resp_no_inf==1),sum(will_imp_resp_no_inf==0), sum(isnan(will_imp_resp_no_inf))])
title('Infusion')

imp_resp_inf = data.ImprovedRespNum(infusion_trials);
imp_resp_no_inf = data.ImprovedRespNum(~infusion_trials);
figure(2)
subplot(1,2,1)
bar([sum(imp_resp_inf==1),sum(imp_resp_inf==0), sum(isnan(imp_resp_inf))])
subplot(1,2,2)
bar([sum(imp_resp_no_inf==1),sum(imp_resp_no_inf==0), sum(isnan(imp_resp_no_inf))])
title('No Infusion')