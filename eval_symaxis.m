nimages = 141;

precision = zeros(1,10);
recall = zeros(1,10);
for ntrials = 1:10
    tp = 0;
    fp = 0;
    fn = 0;
    for im_index = 1:nimages

        load(sprintf('SymDataOurs/Model4/gt%03d.mat',im_index));
        load(sprintf('SymDataOurs/Model4/dt%03d.mat',im_index));

        angle0 = gt{1};
        displ0 = gt{2};
        midpoint0 = gt{3};
        seglen0 = gt{4};

        angle = dt{1};
        displ = dt{2};
        cnfdc = dt{3};

        nhits = 0;
        for i = 1:ntrials
            ddispl = abs(signeddistlinepoint(angle(i),displ(i),midpoint0));
            angleI = angle(i)/pi*180;
            angle0 = angle0/pi*180;
            dangle = min(abs(angleI-angle0),180-abs(angleI-angle0));
            if dangle < 10 && ddispl < 0.2*seglen0
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
plot(recall,precision)
xlabel('recall'), ylabel('precision'), axis([0 1 0 1])
title('single')

% save recall4 recall
% save precision4 precision