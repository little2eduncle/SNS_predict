data_dir = 'D:\Documents\MATLAB\SNS\sim_data';

S_ZZ = load(fullfile(data_dir, 'C_ZZ_intensity_all_10.mat'));
S_XX = load(fullfile(data_dir, 'C_XX_intensity_all_10.mat'));

L_values = S_ZZ.L_values;
num_runs = numel(L_values);

for idx = 1:400
    L = L_values(idx);

    % ---------- ZZ ----------
    C_ZZ_counts3 = S_ZZ.C_ZZ_counts3_all(:, :, :, idx);

    % ---------- XX ----------
    C_XX_counts3   = S_XX.C_XX_counts3_all(:, :, :, idx);

    
    fprintf('L = %.1f km\n', L);
    
    disp('C_XX_counts3 =');
    disp(C_XX_counts3);

    % 如果你后面要在这里继续处理，就直接用上面这些变量
end