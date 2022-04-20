clear, close all

T = readtable('Analysis_cases_allN.xlsx');

countries = {'Austria','Belgium','Bulgaria','Croatia','Czech_Republic',...
    'Estonia','Finland','France','Germany','Greece',... 
    'Hungary','Ireland','Italy','Latvia','Lithuania',... 
    'Netherlands','Poland','Portugal','Romania','Slovakia',...
    'Slovenia','Spain','Sweden','Switzerland','United_Kingdom'};

country = string(T.Country); date = datenum(T.Date); suma7 = table2array(T(1:end,3:24)); suma14 = table2array(T(1:end,25:46)); suma21 = table2array(T(1:end,47:68));
N = [11:17]; ind_N = [5:11]; NORM7 = []; NORM14 = []; NORM21 = []; jump_N_con = zeros(length(country),1); jump_glob = zeros(length(country),1); jump_union = zeros(length(country),1); 
count = 1; color = jet(length(countries)); e14 = table2array(T(1:end,end));

for i = 1:length(countries)
    v_ind = find(country == countries{i});
    DATES = date(v_ind); SUMA7 = suma7(v_ind,ind_N); SUMA14 = suma14(v_ind,ind_N); SUMA21 = suma21(v_ind,ind_N);
    norm7 = SUMA7./SUMA7(:,4); norm14 = SUMA14./SUMA14(:,4); norm21 = SUMA21./SUMA21(:,4);
    NORM7 = [NORM7; norm7]; NORM14 = [NORM14; norm14]; NORM21 = [NORM21; norm21];
    colors = jet(length(DATES)); 
    
    for j = 1:length(DATES)

        casos = SUMA7(j,2:6); max7 = max(casos); min7 = min(casos); diff_all = (max7-min7)/max7; 
        diff7 = abs(diff(casos)./casos(1:4));
        if max(diff7) > 0.25 % consecutive N's jump, around N14
            jump_N_con(count) = 1;
            jump_union(count) = 1;
        end
        if diff_all > 0.35 % global jump, around N14
            jump_glob(count) = 1;
            jump_union(count) = 1;
        end
        
        count = count + 1;
    end
    
    hold off
end

sum(jump_N_con)
sum(jump_glob)
sum(jump_union)


% Write Excel file indicating unstable predictions
filename = 'Analysis_stability.xlsx';
write_instability(filename,country,date,NORM7,NORM14,NORM21,jump_N_con,jump_glob,jump_union,e14)


% Find unstable countries at every date for Patterns_predictions.m 
ind = find(jump_union == 1); DATE_in = date(ind); COUNTRY_in = country(ind); e14_in = e14(ind);


% Figure instability
e14_ins = e14.*jump_union; e14_ins = nonzeros(e14_ins);
Nbins = 14; edges = linspace(-1.3,1.3,Nbins); 
edges_ticks_s = {};
for i = 1:length(edges)-1
    edges_ticks_s{end+1} = [num2str(round(edges(i),3)), ' - ', num2str(round(edges(i+1),3))];
end


f1001 = figure(1001);
hINT106 = histogram(e14_ins,edges);
xlabel('ER14')
ylabel('Prediccions inestables')
ax=hINT106.Parent;
ax.XMinorTick = 'off'; %more tick marks
edges_ticks = edges+(edges(2)-edges(1))/2; edges_ticks = edges_ticks(1:end-1);
set(ax, 'XTick', edges_ticks) 
set(gca, 'XtickLabel', edges_ticks_s);
xtickangle(60)

fs = figure(1000);
subplot(2,2,1)
T = readtable('Errors_N.xlsx'); x = 7:28;
errors = table2array(T(1:end,3:end));
mean_abs = mean(abs(errors)); mean_abs_7 = mean_abs(1:22); mean_abs_14 = mean_abs(23:44); mean_abs_21 = mean_abs(45:66);

plot(x,mean_abs_7,'color',[0 0.7143 1.0000],'linewidth',2)
hold on 
plot(x,mean_abs_14,'color',[0 0.9 0.3],'linewidth',2)
plot(x,mean_abs_21,'color',[1.0000 0.1429 0],'linewidth',2)
xlim([6.5 28.5])
ylim([0 0.5])
legend('7 days','14 days','21 days','location','northeast')
xlabel('Number of days used for fitting')
ylabel('Mean absolute error')
set(gca, 'FontName', 'DejaVu Sans')
hold off

subplot(2,2,2)
T = readtable('Analysis_cases_allN.xlsx'); country = string(T.Country); date = datenum(T.Date); suma7 = table2array(T(1:end,3:24)); suma14 = table2array(T(1:end,25:46)); suma21 = table2array(T(1:end,47:68));

N = 11:17; ind_N = 5:11; 
count = 1; 
v_ind = find(country == 'Slovenia');
DATES = date(v_ind); SUMA7 = suma7(v_ind,ind_N); 
colors = jet(length(DATES)); 
for j = 1:length(DATES)
    plot(N,SUMA7(j,:)/SUMA7(j,4), '.-', 'MarkerSize',20,'color',colors(j,:))
    xlabel('Number of days used for fitting')
    ylabel('Normalized sum of predicted cases (next 7 days)')
    xlim([10.5 17.5])
    hold on
end
set(gca, 'FontName', 'DejaVu Sans')
ylim([0 1.6])
hold off

subplot(2,2,3) 
hINT102 = histogram(e14,edges,'facecolor',[0 20 185]/255);
hold on
histogram(e14_ins,edges,'facealpha',.8,'facecolor',[255 0 0]/255);
xlabel('14-day cumulative relative error')
ylabel('Predictions')
ax=hINT102.Parent;
ax.XMinorTick = 'off'; 
edges_ticks = edges+(edges(2)-edges(1))/2; edges_ticks = edges_ticks(1:end-1);
legend({' All predictions',strcat(" Unstable", string(newline), " predictions")},'location','northwest')
set(ax, 'XTick', edges_ticks) 
set(gca, 'XtickLabel', edges_ticks_s);
set(gca, 'FontName', 'DejaVu Sans')
xtickangle(60)
hold off

subplot(2,2,4)
ratio = (hINT106.Values) ./ hINT102.Values * 100;
bar(edges_ticks,ratio,1,'facecolor',[0 20 185]/255,'FaceAlpha', 0.6);
xlabel('14-day cumulative relative error')
ylabel('Percentage of unstable predictions')
ax.FontSize = 12;
ax=hINT102.Parent;
ax.XMinorTick = 'off'; 
set(gca, 'XTick', edges_ticks) 
set(gca, 'XtickLabel', edges_ticks_s);

xtickangle(60)
set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
set(gca, 'FontName', 'DejaVu Sans')


%%
function [] = write_instability(filename,country,date,NORM_7,NORM_14,NORM_21,jump_N_con,jump_glob,jump_union,e14)

pagina = 'Hoja1';
xlswrite(filename,string(datestr(date)),pagina,'A2')
xlswrite(filename,country,pagina,'B2')
xlswrite(filename,round(NORM_7,3),pagina,'C2')
xlswrite(filename,jump_N_con,pagina,'J2')
xlswrite(filename,jump_glob,pagina,'K2')
xlswrite(filename,jump_union,pagina,'L2')
xlswrite(filename,e14,pagina,'M2')

xlswrite(filename,{'Date'},pagina,'A1')
xlswrite(filename,{'Country'},pagina,'B1')
xlswrite(filename,{'SUMA7 11'},pagina,'C1')
xlswrite(filename,{'SUMA7 12'},pagina,'D1')
xlswrite(filename,{'SUMA7 13'},pagina,'E1')
xlswrite(filename,{'SUMA7 14'},pagina,'F1')
xlswrite(filename,{'SUMA7 15'},pagina,'G1')
xlswrite(filename,{'SUMA7 16'},pagina,'H1')
xlswrite(filename,{'SUMA7 17'},pagina,'I1')
xlswrite(filename,{'Jump consecutive N'},pagina,'J1')
xlswrite(filename,{'Jump global'},pagina,'K1')
xlswrite(filename,{'Any jump'},pagina,'L1')
xlswrite(filename,{'E14'},pagina,'M1')

end
