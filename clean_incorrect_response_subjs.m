function clean_incorrect_response_subjs

%Where are the subjects data located
data_path = 'E:\Users\wilsonj3\Google Drive\NFB_response\IncorrectResponse';


%Grab dir list from datalocation
expression = '.*\d{3,9}.*';
data_dirs = regexp(glob([data_path '/SON*']), expression, 'match'); %Think of a way to just get the uniue ids here , probably with cellfun
data_dirs = [data_dirs{:}];


if ~exist([data_path '/cleaned_already.txt'],'file')
    for data_dir = data_dirs
        files = glob([data_dir{:} '*Run_*.csv'])';
        for file = files
            %read in data
            T = readtable(file{:});
            
            %Replace response strings
            T.WillImpRespText = cellfun(@change_yes_no, T.WillImpRespText, 'UniformOutput', false);
            T.ImprovedRespText = cellfun(@change_yes_no, T.ImprovedRespText, 'UniformOutput', false);
            
            %Resave the files to original place
            writetable(T,file{:})
        end
    end
else
    fprintf('Files have already been cleaned! Exiting...')
end



function y = change_yes_no(str)
if strcmpi('yes',str)
    y = 'No';
elseif strcmpi('no',str)
    y = 'Yes';
else
    y = str;
end


%What were you thinking? The motor inputs are the motor inputs just
%the text (ie responses need to be switched)
% % %         %Save these for original motor regressor creation
% % %         tmp_WillImpReapNum= T.WillImpRespNum;
% % %         tmp_ImprovedRespNum = T.ImprovedRespNum;
% % %
% % %         %Swap the inputs
% % %         T.WillImpRespNum=swap_inputs(T.WillImpRespNum);
% % %         T.ImprovedRespNum=swap_inputs(T.ImprovedRespNum);

% % function data=swap_inputs(data)
% % %Swaps 2 for 7 and vice versa
% %
% % two_idx = data==2;
% % seven_idx = data==7;
% % data(two_idx)=7;
% % data(seven_idx)=2;
