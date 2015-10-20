clearvars
close all
%%

image = 'I17*';

input_path = '~/workspace/SymmetryDBpp/S/';
output_path = '~/workspace/SymmetryDBpp/output/S/';
gt_path = '~/workspace/SymmetryDBpp/S/';

%%
output_files = batch_find_symmetry([input_path image '.png'], output_path)

%%

figure(1)

image = 'I171';
output_path = '~/workspace/SymmetryDBpp/output/S/sigma_2/';

subplot(1,2,1)
show_result([gt_path image '.mat'], [gt_path image '.png']);

subplot(1,2,2)
% show_result(output_files{1});
show_result([output_path image '.mat']);
