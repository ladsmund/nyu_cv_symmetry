clearvars

THRESHOLD_ANGLE = 10 * pi / 180;
THRESHOLD_DISTANCE = @(l1, l2) 0.2 * min(l1,l2);

results_folder = '~/workspace/SymmetryDBpp/output/S/sigmas_02_alpha_10/';
results_folder = '~/workspace/SymmetryDBpp/output_batch/S/';
gt_folder = '~/workspace/SymmetryDBpp/S/';


result_folders = dir('~/workspace/SymmetryDBpp/output/S/');
myDir = find(vertcat(result_folders.isdir));
myDir = myDir(3:end);
%%
% result_folders = result_folders(myDir);

convert_segments = @(s) f_gt(s(1),s(3),s(2),s(4));

%%

figure
hold on
%%
plotNum = 1;
legends = cell(numel(myDir),1);
for dirIndx = myDir'
        
    folder = result_folders(dirIndx).name
    disp([results_folder folder '/*.mat'])
    files = dir([results_folder folder '/*.mat']);
    % files = files(1);

    nimages = numel(files);
    ntrials = 10;

    %%
    precision = zeros(1,10);
    recall = zeros(1,10);

    fail_angle = 0;
    fail_center = 0;

    for ntrials = 1:10
        fprintf('ntrials: %i\n',ntrials);
        tp = 0;
        fp = 0;
        fn = 0;
        for im_index = 1:nimages
            res = load([results_folder folder '/' files(im_index).name],'segments');
            gt  = load([gt_folder files(im_index).name],'segments');        

            [angle0, displ0, midpoin0, seglen0] = convert_segments(gt.segments{1});

            nhits = 0;
            for i = 1:ntrials            
                [angle,displ,midpoint,seglen] = convert_segments(res.segments{i});
                if THRESHOLD_DISTANCE(seglen0,seglen) < norm(midpoin0 - midpoint)
                    fail_center = fail_center + 1;
                elseif THRESHOLD_ANGLE < (max(angle0,angle) - min(angle0,angle))                
                    fail_angle = fail_angle + 1;
                else
                    nhits = nhits+1;
                end
            end
            if nhits > 0;
                tp = tp+1;
            else
                fn = fn+1;
            end
            fp = fp+(ntrials-nhits);
        end
        precision(ntrials) = tp/(tp+fp);
        recall(ntrials) = tp/(tp+fn);
    %     [tp/(tp+fp) tp/(tp+fn)]
    %     return
    end
    plot(recall,precision);
    legends{plotNum} = folder;
    plotNum = plotNum + 1;
end
%%
hold off
xlabel('recall'), ylabel('precision'), axis([0 1 0 1])
title('single')    
legend(legends);
    
    
    
