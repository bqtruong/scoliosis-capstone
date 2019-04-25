function [] = plot_mediastinum()
    close all force
    
    % Load and combine variables

    load('full_stats.mat');
    % load('full_ms.mat');
    load('full_lungs.mat');
    ages_str = a;
    ages = zeros(length(ages_str),1);
    for x = 1:length(ages_str)
        if ages_str{x}(end) == 'M'
            ages(x) = str2double(ages_str{x}(1:3)) / 12;
        elseif ages_str{x}(end) == 'Y'
            ages(x) = str2double(ages_str{x}(1:3));
        else
            disp("Something's wrong!");
        end
    end
    volumes = v;
    gender = g;
    height = h;
    weight = w;
    lv = l;
    
    % Calculate BMI
    BMI = weight ./ (height .^ 2);
    
    % Outlier adjustment
    height(ages == 116) = [];
    weight(ages == 116) = [];
    BMI(ages == 116) = [];
    volumes(ages == 116) = [];
    lv(ages == 116) = [];
    ages(ages == 116) = [];
    height(40) = [];
    weight(40) = [];
    BMI(40) = [];
    volumes(40) = [];
    lv(40) = [];
    ages(40) = [];    
    
    f1 = figure(1);
    scatter(BMI, volumes)
    title('BMI vs. Mediastinal Volume')
    ylabel('Mediastinal Volume (mm^3)')
    xlabel('BMI')
    saveas(f1, 'BMI.png');
    f2 = figure(2);
    scatter(height(height ~= 0), volumes(height ~= 0))
    title('Height vs. Mediastinal Volume')
    ylabel('Mediastinal Volume (mm^3)')
    xlabel('Height (m)')
    saveas(f2, 'Height.png');
    f3 = figure(3);
    scatter(weight(weight ~= 0), volumes(weight ~= 0))
    title('Weight vs. Mediastinal Volume')
    ylabel('Mediastinal Volume (mm^3)')
    xlabel('Weight (kg)')
    saveas(f3, 'Weight.png');
    f4 = figure(4);
    s = scatter(height(height ~= 0 & weight ~= 0), weight(height ~= 0 & weight ~= 0), volumes(height ~= 0 & weight ~= 0)/max(volumes(height ~= 0 & weight ~= 0))*100, volumes(height ~= 0 & weight ~= 0), 'filled');
    sc = colorbar('eastoutside');
    sc.Label.String = 'Mediastinal Volume (mm^3)';
    alpha(s, 0.8);
    xlabel('Height (m)');
    ylabel('Weight (kg)');
    title('Mediastinal Volume (mm^3) Across Height (m) and Weight (kg)')
    saveas(f4, 'Height_Weight_Volume.png');
    f5 = figure(5);
    scatter(ages, volumes)
    title('Age (yr) vs. Mediastinal Volume (mm^3)');
    xlabel('Age (yr)');
    ylabel('Mediastinal Volume (mm^3)');
    saveas(f5, 'Age.png');
    f6 = figure(6);
    scatter3(height(height ~= 0 & weight ~= 0), weight(height ~= 0 & weight ~= 0), ages(height ~= 0 & weight ~= 0), volumes(height ~= 0 & weight ~= 0)/max(volumes(height ~= 0 & weight ~= 0))*100, volumes(height ~= 0 & weight ~= 0));
    sc = colorbar('eastoutside');
    sc.Label.String = 'Mediastinal Volume (mm^3)';
    title('Mediastinal Volume Across Height, Weight, and Age');
    xlabel('Height (m)');
    ylabel('Weight (kg)');
    zlabel('Age (yr)');
    saveas(f6, 'Height_Weight_Age_Volume.png');
    saveas(f6, 'Height_Weight_Age_Volume.fig');
    f7 = figure(7);
    scatter3(height(height ~= 0 & weight ~= 0), weight(height ~= 0 & weight ~= 0), volumes(height ~= 0 & weight ~= 0), lv(height ~= 0 & weight ~= 0)/max(lv(height ~= 0 & weight ~= 0))*100, lv(height ~= 0 & weight ~= 0));
    sc = colorbar('eastoutside');
    sc.Label.String = 'Lung Volume (mm^3)';
    title('Lung Volume Across Height, Weight, and Mediastinal Volume');
    xlabel('Height (m)');
    ylabel('Weight (kg)');
    zlabel('Mediastinal Volume (mm^3)');
    saveas(f7, 'Height_Weight_MS_Lung.png');
    saveas(f7, 'Height_Weight_MS_Lung.fig');
    f8 = figure(8);
    scatter(volumes, lv);
    title('Mediastinal Volume vs. Lung Volume');
    xlabel('Mediastinal Volume (mm^3)');
    ylabel('Lung Volume (mm^3)');
    saveas(f8, 'MS_Lung.png');
    
    f9 = figure(9);
    
    subplot(2,2,1);
    scatter(ages, volumes)
    title('Age (yr) vs. Mediastinal Volume (mm^3)', 'Units', 'normalized', 'Position', [0.5, 1.1, 0]);
    xlabel('Age (yr)');
    ylabel('Mediastinal Volume (mm^3)');
    
    subplot(2,2,2);
    scatter(height(height ~= 0), volumes(height ~= 0))
    title('Height vs. Mediastinal Volume', 'Units', 'normalized', 'Position', [0.5, 1.1, 0])
    ylabel('Mediastinal Volume (mm^3)')
    xlabel('Height (m)')

    subplot(2,2,3);
    scatter(weight(weight ~= 0), volumes(weight ~= 0))
    title('Weight vs. Mediastinal Volume', 'Units', 'normalized', 'Position', [0.5, 1.1, 0])
    ylabel('Mediastinal Volume (mm^3)')
    xlabel('Weight (kg)')
    saveas(f9, 'AllParams.png');
end