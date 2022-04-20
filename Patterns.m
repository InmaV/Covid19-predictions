clear, close all

warning('off','all')

T = readtable('Data_EU+EFTA+UK_2020.xlsx');
[~,~,popul] = xlsread('population.xlsx');

countries = {'Austria','Belgium','Bulgaria','Croatia','Czech_Republic',...
    'Finland','France','Germany','Greece',... 
    'Hungary','Ireland','Italy','Lithuania',... 
    'Netherlands','Poland','Portugal','Romania','Slovakia',...
    'Slovenia','Spain','Sweden','Switzerland','United_Kingdom'};
countries2 = {'Austria','Belgium','Bulgaria','Croatia','Czech Republic',...
    'Finland','France','Germany','Greece',... 
    'Hungary','Ireland','Italy','Lithuania',... 
    'Netherlands','Poland','Portugal','Romania','Slovakia',...
    'Slovenia','Spain','Sweden','Switzerland','United Kingdom'};

pop = zeros(length(countries),1);
for k = 1:length(countries)
    pop(k) = popul{strcmp(popul(:,1),countries(k)),2}*1e3;
end
countries{end+1} = 'Europe'; countries2{end+1} = 'Europe'; pop = [pop; 527862990];

PatternTime = 91;
colors = jet(length(countries));
initial_date = datenum('01/09/2020','dd/mm/yyyy');
final_date = datenum('30/11/2020','dd/mm/yyyy');
times = initial_date:final_date;
T.Day = datenum(T.Day);
alph = 0.025;
sigma = zeros(length(countries),1);
no_report = zeros(length(countries),1);
diff_weight = zeros(length(countries),1);
p_weight = zeros(length(countries),length(times));
W_matrix = zeros(length(countries),7);

for ID = 1:length(countries)
    
    color = colors(ID,:);
    eval(['cas = T.' countries{ID} ';'])
    dia_k = T.Day(2:end-3);
    count = 0;

    dies_p = datenum(dia_k(2:end));
    cum = cas(2:end);
    new_k = cas(2:end)-cas(1:end-1);
    new_k7 = filter(ones(1,7)/7,1,new_k);
    new_k(end-2:end) = [];
    new_k7(1:3) = [];
    POS = zeros(length(dia_k),1);
    POS(1:7) = [1 3 5 2 4 5 1]';

    for k = 8:length(POS)
        if POS(k-7) == 5
            POS(k) = 1;
        else
            POS(k) = POS(k-7)+1;
        end
    end

    f1 = figure(1);
    clf
    dia_ki  = dia_k(dia_k<=final_date); 
    new_ki  = new_k(dia_k<=final_date); 
    new_k7i = new_k7(dia_k<=final_date); 
    W = reporting_pattern(dia_ki,new_ki,PatternTime);
    W_matrix(ID,:) = W;

    for j = 1:length(dia_ki)
        if j > length(dia_ki) - PatternTime
            wd = weekdayC(dia_ki(j));
            col = color;
            if new_ki(j) == 0
                if j == length(dia_ki)
                    valor = new_ki(j-1);
                    i1 = weekdayC(dia_ki(j)); i2 = weekdayC(dia_ki(j-1));
                    datestr(dia_ki(j))
                    new_ki(j) = valor*(W(i1)/(W(i1)+W(i2))); 
                    new_ki(j-1) = valor*(W(i2)/(W(i1)+W(i2)));
                else
                    valor = new_ki(j+1);
                    i1 = weekdayC(dia_ki(j)); i2 = weekdayC(dia_ki(j+1));
                    datestr(dia_ki(j))
                    new_ki(j) = valor*(W(i1)/(W(i1)+W(i2))); 
                    new_ki(j+1) = valor*(W(i2)/(W(i1)+W(i2))); 
                end
            end
            plot(wd+0.1*(POS(j)-3),new_ki(j)/new_k7i(j),'o','Color',col,'MarkerFaceColor',col*0.5+0.5,'LineWidth',1.5)
            hold on
            if new_k7i(j) == 0
                p_weight(ID,j-(length(dia_ki) - PatternTime)) = 0;
            else
                p_weight(ID,j-(length(dia_ki) - PatternTime)) = new_ki(j)/new_k7i(j);
            end
        end
    end

    for k=1:7
        plot(k+[-0.4 0.4],W(k)*[1 1],'-','LineWidth',1.5,'Color',0.8*color)
    end
    axis([0.5 7.5 0 3.5])
    ax=gca;
    ax.XTick = 1:7;
    ax.XTickLabel = {'Mo.','Tu.','We.','Th.','Fr.','Sa.','Su.'};
    ax.LineWidth = 1.3;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ylabel('Weight')
    title(countries2{ID})
    
    weights_country = p_weight(ID,:);
    w_monday = weights_country(weekdayC(times)==1);
    w_tuesday = weights_country(weekdayC(times)==2);
    w_wednesday = weights_country(weekdayC(times)==3);
    w_thursday = weights_country(weekdayC(times)==4);
    w_friday = weights_country(weekdayC(times)==5);
    w_saturday = weights_country(weekdayC(times)==6);
    w_sunday = weights_country(weekdayC(times)==7);
    
    sigma(ID) = mean([std(w_monday) std(w_tuesday) std(w_wednesday) std(w_thursday) std(w_friday) std(w_saturday) std(w_sunday)]);
    diff_weight(ID) = max(W) - min(W);
    no_report(ID) = count;    
end

lpop = log10(pop);

% Write Excel
filename = 'Pattern_analysis.csv';
writecsv(filename,countries',lpop,sigma,diff_weight)

% % Remove countries with more than 2 days without report
% dies_no_report = 3; v = find(no_report < dies_no_report); 
% for i = 1:length(v)
%     countries_restants{i} = countries{v(i)};
% end
% pop = pop(no_report < dies_no_report);
% sigma = sigma(no_report < dies_no_report);
% diff_weight = diff_weight(no_report < dies_no_report);
% 
% % Remove countries with less than 1 million inhabitants:
% lim_pop = 1e6;
% v = find(pop > lim_pop);
% for i = 1:length(v)
%     countries_r{i} = countries_restants{v(i)};
% end
% sigma = sigma(pop > lim_pop);
% diff_weight = diff_weight(pop > lim_pop);
% pop = pop(pop > lim_pop); 
% 
% % Plots 
% f2 = figure(2);
% scatter3(pop,sigma,diff_weight,'filled')
% xlabel('Population'); ylabel('Standard deviation'); zlabel('Difference max-min')
% grid on
% 
% f3 = figure(3);
% scatter3(log(pop),sigma,diff_weight,'filled')
% xlabel('log(Population)'); ylabel('Standard deviation'); zlabel('Difference max-min')
% grid on
% 
% f4 = figure(4);
% scatter3((pop-mean(pop))./std(pop),(sigma-mean(sigma))./std(sigma),(diff_weight-mean(diff_weight))./std(diff_weight),'filled')
% xlabel('Population'); ylabel('Standard deviation'); zlabel('Difference max-min')
% grid on
% 
% f5 = figure(5);
% scatter3((log(pop)-mean(log(pop)))./std(log(pop)),(sigma-mean(sigma))./std(sigma),(diff_weight-mean(diff_weight))./std(diff_weight),'filled')
% xlabel('Population'); ylabel('Standard deviation'); zlabel('Difference max-min')
% grid on

%%
function [] = writecsv(filename,countries,pop,sigma,diff_weight)

fid = fopen(filename, 'w+');
for iD = 1:length(countries)
    if iD < length(countries)
        fprintf(fid, '%s,', countries{iD});
    else
        fprintf(fid, '%s\n', countries{iD});
    end
end

for iD = 1:length(countries)
    if iD < length(countries)
        fprintf(fid, '%d,', round(pop(iD),4));
    else
       fprintf(fid, '%d\n', round(pop(iD),4));
    end
end

for iD = 1:length(countries)
    if iD < length(countries)
        fprintf(fid, '%d,', round(sigma(iD),4));
    else
       fprintf(fid, '%d\n', round(sigma(iD),4));
    end
end

for iD = 1:length(countries)
    if iD < length(countries)
        fprintf(fid, '%d,', round(diff_weight(iD),4));
    else
       fprintf(fid, '%d\n', round(diff_weight(iD),4));
    end
end

fclose(fid);

end
function [w,C] = reporting_pattern(dies,casos,PatternTime)

casosf = filter(ones(1,7)/7,1,casos);
diesf = dies - 3;

diesf = diesf(end-PatternTime+1:end);
casosf = casosf(end-PatternTime+1:end);

C = cell(1,7);
for k = 1:PatternTime
    if casosf(k)>0
        d = weekday(diesf(k))-1;
        d(d==0) = 7;
        C{d} = [ C{d} casos(dies==diesf(k))/casosf(k)];
    end
end

w = ones(7,1);
for k = 1:7
    if ~isempty(C{k})
        w(k) = mean(C{k});
    end
end

end
function wd = weekdayC(dia)
wd = weekday(dia)-1;
wd(wd==0) = 7;
end
