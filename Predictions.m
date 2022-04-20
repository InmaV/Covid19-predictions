% clear, close all

warning('off','all')

T = readtable('Data_EU+EFTA+UK_2020.xlsx');
[~,~,popul] = xlsread('population.xlsx');

countries = {'Austria','Belgium','Bulgaria','Croatia','Czech_Republic',...
    'Finland','France','Germany','Greece',... 
    'Hungary','Ireland','Italy','Lithuania',... 
    'Netherlands','Poland','Portugal','Romania','Slovakia',...
    'Slovenia','Spain','Sweden','Switzerland','United_Kingdom'};

pop = zeros(length(countries),1);
for k = 1:length(countries)
    pop(k) = popul{strcmp(popul(:,1),countries(k)),2}*1e3;
end

PatternTime = 91;
initial_date = datenum('01/09/2020','dd/mm/yyyy');
final_date = datenum('30/11/2020','dd/mm/yyyy');
times1 = datenum('05/09/2020','dd/mm/yyyy'):7:final_date;
times2 = datenum('01/09/2020','dd/mm/yyyy'):7:final_date;
times = union(times1,times2);
T.Day = datenum(T.Day);
alph = 0.025;

constant_v = [true true true true];
weights_v  = [false false true true];
newfit_v   = [false true true false];
colors = [[0 0.62 0];[0 0 1];[1 0 0];[1 0.75 0]];

Ndays_predicted = 21;
weightt = 1;

NAME = {'B','F','H','I'};
country_error = length(countries);

all_days = initial_date:final_date; v_holidays = holidays_countries(all_days,countries);
  
v_lim = linspace(0.2,0.6,21); v_error_total = []; v_error_holidays = []; v_error_no_holidays = [];
count_my_100 = 0; count_holidays = 0;
fest = [];
v_Nfitting = [14]; % linspace(7,28,22); %linspace(12,16,5);

ERRORS_7 = zeros(length(v_Nfitting),length(countries)*length(times));
ERRORS_14_ac = zeros(length(v_Nfitting),length(countries)*length(times)); 
ERRORS_21_ac = zeros(length(v_Nfitting),length(countries)*length(times));
SUM_CASOS_21 = zeros(length(v_Nfitting),length(countries)*length(times)); 
SUM_CASOS_14 = zeros(length(v_Nfitting),length(countries)*length(times)); 
SUM_CASOS_7= zeros(length(v_Nfitting),length(countries)*length(times)); 
COUNTRIES = cell([length(countries)*length(times),1]);
DATES_PRED2 = zeros(length(countries)*length(times),1);
v_ins = zeros(length(countries)*length(times),1);


for Nfit = 1:length(v_Nfitting)
count = 1;

for ijk = 3%1:length(weights_v)
    CASOS_TOTALS = zeros(length(countries),Ndays_predicted,length(times));
    
    weights  = weights_v(ijk);
    constant = constant_v(ijk);
    newfit   = newfit_v(ijk);
    color = colors(ijk,:);

    if newfit
        F_value = 1;
    else
        F_value = 0;
    end
    
    MAT_r7a = zeros(length(times),length(countries)); MAT_r7my5 = zeros(length(times),length(countries));
    MAT_r7my10 = zeros(length(times),length(countries)); MAT_r7my15 = zeros(length(times),length(countries));
    MAT_r75 = zeros(length(times),length(countries)); MAT_r710 = zeros(length(times),length(countries));
    MAT_r73 = zeros(length(times),length(countries)); MAT_r74 = zeros(length(times),length(countries)); 
    MAT_r7my3 = zeros(length(times),length(countries)); MAT_pa14my3 = zeros(length(times),length(countries));
    MAT_pa14a = zeros(length(times),length(countries)); MAT_pa14my5 = zeros(length(times),length(countries));
    MAT_pa14my10 = zeros(length(times),length(countries)); MAT_pa14my15 = zeros(length(times),length(countries));
    MAT_pa145 = zeros(length(times),length(countries)); MAT_pa1410 = zeros(length(times),length(countries));
    MAT_e7 = zeros(length(times),length(countries)); MAT_e14 = zeros(length(times),length(countries));
    MAT_e21 = zeros(length(times),length(countries)); MAT_e7_r = zeros(length(times),length(countries));
    MAT_e14_r = zeros(length(times),length(countries)); MAT_e21_r = zeros(length(times),length(countries));
    MAT_m = zeros(length(times),length(countries)); MAT_e7_r_ac = zeros(length(times),length(countries));
    MAT_e14_r_ac = zeros(length(times),length(countries)); MAT_e21_r_ac = zeros(length(times),length(countries));
    MAT_r7my6 = zeros(length(times),length(countries)); MAT_pa14my6 = zeros(length(times),length(countries));
    countries_R2 = {}; dia_R2 = {}; setmana_R2 = {}; pa14_R2 = []; diff_R2 = []; rho_R2_6 = []; casos_R2 = []; er7_R2 = [];
    IA14_R2 = []; IA14_R2_tot = []; pa14_R2_tot = []; rho_R2_6_tot = []; casos_R2_tot = []; rho_R2_15 = []; rho_R2_10 = []; rho_R2_4 = [];
    countries_drh = {}; dia_drh = {}; setmana_drh = {}; pa14_drh = []; diff_drh = []; rho_drh = []; casos_drh = []; er7_drh = [];
    IA14_drh = []; countries_R2_tot = {}; dia_R2_tot = {}; setmana_R2_tot = {}; diff_R2_tot = []; rho_R2_15_tot = []; rho_R2_10_tot = []; 
    rho_R2_4_tot = []; er7_R2_tot = []; eliminats = zeros(length(times),length(countries));
    countries_f = {}; dia_f = {}; setmana_f = {}; pa14_f = []; diff_f = []; rho_f = []; casos_f = []; er7_f = []; IA14_f = [];
    countries_R2_4 = {}; dia_R2_4 = {}; setmana_R2_4 = {}; er7_R2_4 = []; IA14_R2_4 = []; diff_R2_4 = []; rho_R2_6_4 = []; rho_R2_15_4 = []; rho_R2_10_4 = []; rho_R2_4_4 = []; casos_R2_4 = [];
    countries_R2_6 = {}; dia_R2_6 = {}; setmana_R2_6 = {}; er7_R2_6 = []; IA14_R2_6 = []; diff_R2_6 = []; rho_R2_6_6 = []; rho_R2_15_6 = []; rho_R2_10_6 = []; rho_R2_4_6 = []; casos_R2_6 = [];

    for ID = 1:length(countries)
        
        v_fest = v_holidays(ID,:);
        e7 = zeros(length(times),1); e14 = zeros(length(times),1); e21 = zeros(length(times),1); v_m = zeros(length(times),1); 
        a14 = zeros(length(times),1); e7_r = zeros(length(times),1); e14_r = zeros(length(times),1); e21_r = zeros(length(times),1);
        e7_r_ac = zeros(length(times),1); e14_r_ac = zeros(length(times),1); e21_r_ac = zeros(length(times),1);
        
        eval(['cas = T.' countries{ID} ';']);
        dia_k = T.Day(2:end);

        dies_p = datenum(dia_k);
        cum = cas(2:end);
        new_k = cas(2:end)-cas(1:end-1);
        new_k7 = filter(ones(1,7)/7,1,new_k);
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
  
        dia_ki  = dia_k(dia_k<=final_date); 
        new_ki  = new_k(dia_k<=final_date);
        
        if weights
            W = reporting_pattern(dia_ki,new_ki,PatternTime);
        else
            W = ones(7,1);
        end

        Wi = W(weekday(dia_ki-1)); 
        casosp = new_ki./Wi;
        casosp(new_ki==0) = 0;
        casosp(Wi==0) = 0;
        ccp_k = cumsum(casosp); 
        CASOS = zeros(Ndays_predicted,length(times));
        
        Wi_tot = W(weekday(dia_k-1)); % with all days
        casosp_tot = new_k./Wi_tot;
        casosp_tot(new_k==0) = 0;
        casosp_tot(Wi_tot==0) = 0;
        ccp_k_tot = cumsum(casosp_tot);
        
        index1 = find(dia_k==times1(1)):7:find(dia_k==times1(end));
        index2 = find(dia_k==times2(1)):7:find(dia_k==times2(end));
        index = union(index1,index2); 
        a14(1:14) = ccp_k_tot(1:14); a14(15:length(ccp_k_tot)) = ccp_k_tot(15:end) - ccp_k_tot(1:end-14);
        mitjana_casos = mitjana(2,casosp_tot);
        c7 = casosp_tot; 
        id = 5:length(casosp_tot)-3;
        rh_pon = 0.*casosp_tot;
        rh_pon(id) = (c7(id+1)+c7(id+2)+c7(id+3))./(c7(id-2)+c7(id-3)+c7(id-4));
        id = 7:length(casosp_tot)-1;
        
        rh7_c = mitjana(2,rh_pon); rh7_c_a = rh7_c(index); MAT_r7a(:,ID) = rh7_c_a; % actual
        rh7_c_my3 = rh7_c(index-3); MAT_r7my3(:,ID) = rh7_c_my3; % 3 days ago
        rh7_c_my5 = rh7_c(index-5); MAT_r7my5(:,ID) = rh7_c_my5; % 5 days ago
        rh7_c_my6 = rh7_c(index-6); MAT_r7my6(:,ID) = rh7_c_my6; % 6 days ago
        rh7_c_3 = rh7_c(index+3); MAT_r73(:,ID) = rh7_c_3; % 3 days forward
        rh7_c_4 = rh7_c(index+4); MAT_r74(:,ID) = rh7_c_4; % 4 days forward
        rh7_c_5 = rh7_c(index+5); MAT_r75(:,ID) = rh7_c_5; % 5 days forward
        rh7_c_10 = rh7_c(index+10); MAT_r710(:,ID) = rh7_c_10; % 10 days forward
        rh7_c_my15 = rh7_c(index-15); MAT_r7my15(:,ID) = rh7_c_my15; % 15 days ago
        rh7_c_my10 = rh7_c(index-10); MAT_r7my10(:,ID) = rh7_c_my10; % 10 days ago
        diffrh7 = rh7_c_my6 - rh7_c_my15;       
        
        dies_p2 = datenum(dia_k);
        Wi_plot = W(weekday(dies_p2-1));
        casos_plot = new_k./Wi_plot;
        casos_plot(new_k==0) = 0;
        casos_plot(Wi_plot==0) = 0;

        for i_dia = 1:length(times)
            FEST = 0; CEN = 0;
            Nfitting = v_Nfitting(Nfit); 

            xx = dia_ki(find(dia_ki==times(i_dia))-Nfitting+1:find(dia_ki==times(i_dia))); 
            yy = ccp_k(find(dia_ki==times(i_dia))-Nfitting+1:find(dia_ki==times(i_dia))); 
            v_c = casosp(find(dia_ki==times(i_dia))-20+1:find(dia_ki==times(i_dia)));
            
            index = find(dia_ki<=times(i_dia) & dia_ki>=initial_date);
            AC_PATTERN = cum(index); DATES = dia_ki(index); NOUS_PATTERN = new_ki(index); NOUS_PON = casosp(index); AC_PON = ccp_k(index);
            
            Nt = 1e5;
            rng(12345); % seed 
            lim_inf = pop(ID)/2030; 
            w = ones(Nfitting,1);
            [p,m] = fitGompertz(yy,w,Nt,pop(ID),constant,F_value,lim_inf);

            t = xx(end) + (0:Ndays_predicted)'; 
            [g] = gompertz(yy(1),p(1),p(2),p(3),t-xx(1));
            nw = g(2:end)-g(1:end-1); NOUS_PON_PRED = nw; 
            ac = yy(end) + cumsum(nw); AC_PON_PRED = g(2:end);
            t = t(2:end); nw_p = nw; DATES_PRED = t; 
            nw = nw.*W(weekday(t-1)); NOUS_PATTERN_PRED = nw; AC_PATTERN_PRED = yy(end) + cumsum(nw);
            
            t2 = initial_date-1:times(i_dia)+21;
            [g2] = gompertz(yy(1),p(1),p(2),p(3),t2-xx(1));
            nw2 = g2(2:end)-g2(1:end-1); nw2 = nw2'; NOUS_PON_FITTING = nw2;
            ac2 = g2(1) + cumsum(nw2); AC_PON_FITTING = ac2;
            t2 = t2(2:end); nw2 = nw2.*W(weekday(t2-1)); 
            NOUS_PATTERN_FITTING = nw2; 
            AC_PATTERN_FITTING = g2(1) + cumsum(nw2);

            i_d = find(dia_k==times(i_dia)+1);
            e7_r_ac(i_dia) = -((sum(nw(1:7)) - sum(new_k(i_d:i_d+6)))/sum(new_k(i_d:i_d+6))); v_error_total = [v_error_total e7_r_ac(i_dia)]; 
            e14_r_ac(i_dia) = -((sum(nw(1:14)) - sum(new_k(i_d:i_d+13)))/sum(new_k(i_d:i_d+13))); 
            e21_r_ac(i_dia) = -((sum(nw(1:21)) - sum(new_k(i_d:i_d+20)))/sum(new_k(i_d:i_d+20)));

            if mean(new_k(i_d-7:i_d-1)) < 100
                nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                count_my_100 = count_my_100 + 1; eliminats(i_dia,ID) = 1; CEN = 1;
            end
            
            ind_fest = find(all_days == times(i_dia));
            if ind_fest == 1
                if sum(v_fest(ind_fest-0:ind_fest)) > 0 
                    nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                    v_error_holidays = [v_error_holidays e7_r_ac(i_dia)];
                else
                    v_error_no_holidays = [v_error_no_holidays e7_r_ac(i_dia)];
                end
            elseif ind_fest == 2
                if sum(v_fest(ind_fest-1:ind_fest)) > 0
                    nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                    v_error_holidays = [v_error_holidays e7_r_ac(i_dia)];
                else
                    v_error_no_holidays = [v_error_no_holidays e7_r_ac(i_dia)];
                end
            elseif ind_fest == 3
                if sum(v_fest(ind_fest-2:ind_fest)) > 0
                    nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0;
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                    v_error_holidays = [v_error_holidays e7_r_ac(i_dia)];
                else
                    v_error_no_holidays = [v_error_no_holidays e7_r_ac(i_dia)];
                end
            elseif ind_fest == 4
                if sum(v_fest(ind_fest-3:ind_fest)) > 0
                     nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                    v_error_holidays = [v_error_holidays e7_r_ac(i_dia)];
                else
                    v_error_no_holidays = [v_error_no_holidays e7_r_ac(i_dia)];
                end
            elseif ind_fest == 5
                if sum(v_fest(ind_fest-4)) > 0
                    nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                end
            else
                if sum(v_fest(ind_fest-5:ind_fest)) > 0
                    nw(1:end) = 0; e7_r_ac(i_dia) = 0; e14_r_ac(i_dia) = 0; e21_r_ac(i_dia) = 0; 
                    MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
                    count_holidays = count_holidays + 1; eliminats(i_dia,ID) = 1;
                    fest = [fest times(i_dia)]; FEST = 1;
                    v_error_holidays = [v_error_holidays e7_r_ac(i_dia)];
                else
                    v_error_no_holidays = [v_error_no_holidays e7_r_ac(i_dia)];
                end
            end

% %%%%%%% Filter out unstable days %%%%%%%%            
%             for ind_kk = 1:length(COUNTRY_in)
%                 if strcmp(countries{ID},COUNTRY_in{ind_kk}) && times(i_dia) == datenum(DATE_in(ind_kk))
%                     nw(1:end) = 0; e7_r_ac(i_dia) = 0; 
%                     MAT_r7my5(i_dia,ID) = 0; MAT_r75(i_dia,ID) = 0; MAT_r710(i_dia,ID) = 0;
%                     v_ins(count) = 1;
%                 end
%             end

            CASOS(:,i_dia) = nw; SUM_CASOS_7(Nfit,count) = sum(nw(1:7)); SUM_CASOS_14(Nfit,count) = sum(nw(1:14));
            SUM_CASOS_21(Nfit,count) = sum(nw);  ERRORS_7(Nfit,count) = e7_r_ac(i_dia); ERRORS_14_ac(Nfit,count) = e14_r_ac(i_dia); ERRORS_21_ac(Nfit,count) = e21_r_ac(i_dia);

            if Nfit == 1
                COUNTRIES{count} = countries{ID}; DATES_PRED2(count) = times(i_dia);
            end
            count = count + 1;
        end
        
        
        CASOS_TOTALS(ID,:,:) = CASOS;
        MAT_e7_r_ac(:,ID) = e7_r_ac; MAT_e14_r_ac(:,ID) = e14_r_ac; MAT_e21_r_ac(:,ID) = e21_r_ac;
        
    end

% % Write all predictions in a csv file
%     for i = 1:length(times)
%         filename = [datestr(times(i),'yyyymmdd')  '_' NAME{ijk} '.csv'];
%         writecsv(filename,times(i) + (1:Ndays_predicted)',countries,CASOS_TOTALS(:,:,i))
%     end


end

end


v = SUM_CASOS_7(1,:);
index = find(v == 0);
ERRORS_7(:,index) = [];
ERRORS_14_ac(:,index) = [];
ERRORS_21_ac(:,index) = [];
SUM_CASOS_7(:,index) = [];
SUM_CASOS_14(:,index) = [];
SUM_CASOS_21(:,index) = [];
COUNTRIES(index) = [];
DATES_PRED2(index) = [];

% % Write cases and errors in Excel files for different values of Nfitting
% write_analysis_sum('Analysis_cases_allN.xlsx',DATES_PRED2,COUNTRIES,round(SUM_CASOS_7'),round(SUM_CASOS_14'),round(SUM_CASOS_21'),ERRORS_14_ac(8,:)')
% write_analysis_errors('Errors_N.xlsx',DATES_PRED2,COUNTRIES,ERRORS_7',ERRORS_14_ac',ERRORS_21_ac')

% Success rate (in the case only one model and Nfitting are being considered)
['Model ' NAME{ijk}]
accepted_error = [.1 .2 .3 .4 .5]
e_7 = [length(ERRORS_7(abs(ERRORS_7)<=0.1)) length(ERRORS_7(abs(ERRORS_7)<=0.2)) length(ERRORS_7(abs(ERRORS_7)<=0.3)) length(ERRORS_7(abs(ERRORS_7)<=0.4)) length(ERRORS_7(abs(ERRORS_7)<=0.5))]/ length(ERRORS_7)
e_14 = [length(ERRORS_14_ac(abs(ERRORS_14_ac)<=0.1)) length(ERRORS_14_ac(abs(ERRORS_14_ac)<=0.2)) length(ERRORS_14_ac(abs(ERRORS_14_ac)<=0.3)) length(ERRORS_14_ac(abs(ERRORS_14_ac)<=0.4)) length(ERRORS_14_ac(abs(ERRORS_14_ac)<=0.5))]/ length(ERRORS_14_ac)
e_21 = [length(ERRORS_21_ac(abs(ERRORS_21_ac)<=0.1)) length(ERRORS_21_ac(abs(ERRORS_21_ac)<=0.2)) length(ERRORS_21_ac(abs(ERRORS_21_ac)<=0.3)) length(ERRORS_21_ac(abs(ERRORS_21_ac)<=0.4)) length(ERRORS_21_ac(abs(ERRORS_21_ac)<=0.5))]/ length(ERRORS_21_ac)

% Correlation figure (in the case only one model and Nfitting are being considered)
MAT_r7_4_6 = MAT_r74 - MAT_r7my6; 
f27 = figure(27);
vx = []; vy = []; countries_corr = {}; r7_4_6_corr = []; e7_corr = [];
for i = 1:length(countries) 
    e7 = MAT_e7_r_ac(:,i); e7 = e7(e7~=0); r_4 = MAT_r7_4_6(:,i); r_4 = r_4(e7~=0); vy = [vy; e7]; vx = [vx; r_4];
    for j = 1:length(e7)
        countries_corr{end+1} = countries{i};
        r7_4_6_corr = [r7_4_6_corr; r_4(j)]; 
        e7_corr = [e7_corr; e7(j)];
    end
    plot(r_4,e7,'o','MarkerFaceColor',[50 80 225]/255,'MarkerEdgeColor','k') 
    hold on
end
[fitresult, gof] = createFit(vx,vy);
gof = cell2mat(struct2cell(gof)); R2_4_6 = gof(2);
x = linspace(-5,5);
xlabel('Difference in rho')
ylabel('7-day cumulative relative error')
xlim([-2.2 2.2])
ylim([-0.9 0.9])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca, 'FontName', 'DejaVu Sans')
hold off


%%
function [v_holidays] = holidays_countries(dies,countries)
v_holidays = zeros(length(countries),length(dies));
holidays = [datenum('24/09/2020','dd/mm/yyyy') datenum('10/10/2020','dd/mm/yyyy') datenum('26/10/2020','dd/mm/yyyy')...
    datenum('01/11/2020','dd/mm/yyyy') datenum('02/11/2020','dd/mm/yyyy') datenum('11/11/2020','dd/mm/yyyy') ...
    datenum('15/11/2020','dd/mm/yyyy') datenum('22/11/2020','dd/mm/yyyy') datenum('29/11/2020','dd/mm/yyyy')];
v_holidays(1,ismember(dies,holidays)) = 1;

holidays = [datenum('27/09/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy') datenum('11/11/2020','dd/mm/yyyy') ... 
    datenum('15/11/2020','dd/mm/yyyy') datenum('29/11/2020','dd/mm/yyyy')];
v_holidays(2,ismember(dies,holidays)) = 1;

holidays = [datenum('06/09/2020','dd/mm/yyyy') datenum('22/09/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(3,ismember(dies,holidays)) = 1;

holidays = [datenum('08/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(4,ismember(dies,holidays)) = 1;

holidays = [datenum('28/09/2020','dd/mm/yyyy') datenum('28/10/2020','dd/mm/yyyy') datenum('11/11/2020','dd/mm/yyyy') ...
    datenum('17/11/2020','dd/mm/yyyy')];
v_holidays(5,ismember(dies,holidays)) = 1;

holidays = [datenum('10/10/2020','dd/mm/yyyy') datenum('24/10/2020','dd/mm/yyyy') datenum('31/10/2020','dd/mm/yyyy') ...
    datenum('06/11/2020','dd/mm/yyyy') datenum('08/11/2020','dd/mm/yyyy')];
v_holidays(6,ismember(dies,holidays)) = 1;

holidays = [datenum('01/11/2020','dd/mm/yyyy') datenum('11/11/2020','dd/mm/yyyy') datenum('29/11/2020','dd/mm/yyyy')];
v_holidays(7,ismember(dies,holidays)) = 1;

holidays = [datenum('20/09/2020','dd/mm/yyyy') datenum('03/10/2020','dd/mm/yyyy') datenum('04/10/2020','dd/mm/yyyy') ...
    datenum('01/11/2020','dd/mm/yyyy') datenum('15/11/2020','dd/mm/yyyy') datenum('18/11/2020','dd/mm/yyyy') ...
    datenum('22/11/2020','dd/mm/yyyy') datenum('29/11/2020','dd/mm/yyyy')];
v_holidays(8,ismember(dies,holidays)) = 1;

holidays = [datenum('28/10/2020','dd/mm/yyyy')];
v_holidays(9,ismember(dies,holidays)) = 1;

holidays = [datenum('06/10/2020','dd/mm/yyyy') datenum('23/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(10,ismember(dies,holidays)) = 1;

holidays = [datenum('26/10/2020','dd/mm/yyyy') datenum('31/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy') ...
    datenum('02/11/2020','dd/mm/yyyy')];
v_holidays(11,ismember(dies,holidays)) = 1;

holidays = [datenum('19/09/2020','dd/mm/yyyy') datenum('04/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy') ...
    datenum('03/11/2020','dd/mm/yyyy')];
v_holidays(12,ismember(dies,holidays)) = 1;

holidays = [datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(13,ismember(dies,holidays)) = 1;

holidays = [datenum('29/11/2020','dd/mm/yyyy')]; % Vacances tardor 17-25 octubre
v_holidays(14,ismember(dies,holidays)) = 1;

holidays = [datenum('01/11/2020','dd/mm/yyyy') datenum('11/11/2020','dd/mm/yyyy')];
v_holidays(15,ismember(dies,holidays)) = 1;

holidays = [datenum('05/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(16,ismember(dies,holidays)) = 1;

holidays = [datenum('30/11/2020','dd/mm/yyyy')];
v_holidays(17,ismember(dies,holidays)) = 1;

holidays = [datenum('01/09/2020','dd/mm/yyyy') datenum('15/09/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy') ...
    datenum('17/11/2020','dd/mm/yyyy')];
v_holidays(18,ismember(dies,holidays)) = 1;

holidays = [datenum('15/09/2020','dd/mm/yyyy') datenum('31/10/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy') ...
     datenum('11/11/2020','dd/mm/yyyy') datenum('23/11/2020','dd/mm/yyyy')];
v_holidays(19,ismember(dies,holidays)) = 1;

holidays = [datenum('02/09/2020','dd/mm/yyyy') datenum('08/09/2020','dd/mm/yyyy') datenum('11/09/2020','dd/mm/yyyy') ...
    datenum('09/10/2020','dd/mm/yyyy') datenum('12/10/2020','dd/mm/yyyy') datenum('31/10/2020','dd/mm/yyyy') ...
    datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(20,ismember(dies,holidays)) = 1;

holidays = [datenum('24/10/2020','dd/mm/yyyy') datenum('31/10/2020','dd/mm/yyyy') datenum('06/11/2020','dd/mm/yyyy')];
v_holidays(21,ismember(dies,holidays)) = 1;

holidays = [datenum('20/09/2020','dd/mm/yyyy') datenum('01/11/2020','dd/mm/yyyy')];
v_holidays(22,ismember(dies,holidays)) = 1;
end
function [fitresult, gof] = createFit(x,y)

%  Data for fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
% Fit.

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
end
function [] = writecsv(filename,t,countries,cases)

v = cases(:,1);
index = find(v == 0);
cases(index,:) = [];
countries(index) = [];


timeValue = t ; % Time in Matlab format
for iTime = 1:length(timeValue)
   Time{iTime} = datestr(timeValue(iTime), 'yyyy-mm-dd') ; 
end

fid = fopen(filename, 'w+') ;
fprintf(fid, '%s,', 'Dates') ;
for iD = 1:length(countries)
    if iD < length(countries)
        fprintf(fid, '%s,', countries{iD});
    else
        fprintf(fid, '%s\n', countries{iD});
    end
end

for iLine = 1:length(timeValue) % Loop through each time/value row
   fprintf(fid, '%s,', Time{iLine}) ; % Print the time string
   for times = 1:length(countries)
       if times < length(countries)
           fprintf(fid, '%d,', round(cases(times,iLine))); % Print the data values
       else
           fprintf(fid, '%d\n', round(cases(times,iLine))); % Print the data values
       end
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
function [error_total,error_cum,error_new] = ErrorGompertz(cases,K,a,c,w,Fv)

approx = gompertz(cases(1),K,a,c,(0:(length(cases)-1))');
error_cum = sum(abs((cases-approx).*w./cases));

new_cases = cases(2:end)-cases(1:end-1);
new_approx = approx(2:end)-approx(1:end-1);
id = new_cases>0;
new_cases = new_cases(id);
new_approx = new_approx(id);
w = w(2:end); w = w(id);
error_new = 0 + sum(abs((new_cases-new_approx).*w./new_cases));

error_total = error_cum + Fv*error_new;

end
function [g] = gompertz(N0,K,a,c,t)
Nn = N0-c;
g = K*exp(-log(K/Nn)*exp(-a*t))+c;
end
function [p,m] = fitGompertz(yy,w,Nt,pop,ct,Fv,lim_inf)

if ct
    clim = yy(1)-1;
else
    clim = 0;
end

K = [0 pop/5];
A = [0.01 0.5];
C = [lim_inf clim];
LIM = [K' A' C'];
DX = (LIM(2,:)-LIM(1,:))/Nt;
Np = length(DX);
LHS = zeros(Nt,Np);
for i = 1:Nt
    for j = 1:Np
        LHS(i,j) = LIM(1,j) + (i-1+rand)*DX(j);
    end
end
for j=1:Np
    LHS(:,j)=LHS(randperm(Nt),j);
end

id = LHS(:,1)+LHS(:,3) < yy(end) | LHS(:,1)+LHS(:,3) > pop/5;
LHS(id,:) = [];
Nt = size(LHS,1);

et = zeros(Nt,1);
for k = 1:Nt
    et(k) = ErrorGompertz(yy,LHS(k,1),LHS(k,2),LHS(k,3),w,Fv);
end

[m,i] = min(et);
b_guess = LHS(i,:)';

options = optimset('Display','off');
fun = @(b) ErrorGompertz(yy,b(1),b(2),b(3),w,Fv);
p = fmincon(fun, b_guess, [1 0 1; -1 0 -1],[pop;-yy(end)],[],[],[0;0.01;0],[pop/5;0.5;clim],[],options);
m = fun(p);

end
function [w] = weights_function(n,type,dif,w1,l1,w2,l2)

switch type
    case 1
        w = ones(n,1);
    case 2
        w = ones(n,1);
        w(end) = 100;
    case 3
        w = ones(n,1);
        w(end-2:end) = 100;
    case 4
        w = (1:n)';
    case 5
        w = (1:n)';
        w = w.^2;
    case 6
        w = ones(n,1);
        w(1:round(n/2)) = w1;
    case 7
        w = ones(n,1);
        w(round(n/2):end) = w1;
    case 8
        coef = (w1 - 1)/l1*dif + 1;
%         a = ((w1-1)*l2 - (w2-1)*l1)/(l1^2*l2 - l2^2*l1); b = ((w1-1) - a*l1^2)/l1; c = 1;
%         coef = a*dif^2 + b*dif + c;
        w = linspace(coef,1,n)'
%         x = linspace(1,n,n)';
%         c1 = coef; n1 = n; c2 = coef*0.5; n2 = n*0.35;
%         A = ((c2-1)*(n1-1)-(n2-1)*(c1-1))/((n2^2-1)*(n1-1)-(n1^2-1)*(n2-1)); B = (c1-1-(n1^2-1)*A)/(n1-1); C = 1 - A - B;
%         w = A*x.^2 + B*x + C;
    case 9
        coef = (w1 - 1)/l1*dif + 1;
%         a = ((w1-1)*l2 - (w2-1)*l1)/(l1^2*l2 - l2^2*l1); b = ((w1-1) - a*l1^2)/l1; c = 1;
%         coef = a*dif^2 + b*dif + c;
        w = linspace(coef,1,n)'
%         x = linspace(n,1,n)';
%         c1 = coef; n1 = n; c2 = coef*0.5; n2 = n*0.35;
%         A = ((c2-1)*(n1-1)-(n2-1)*(c1-1))/((n2^2-1)*(n1-1)-(n1^2-1)*(n2-1)); B = (c1-1-(n1^2-1)*A)/(n1-1); C = 1 - A - B;
%         w = A*x.^2 + B*x + C;
        
end

end
function [mit_x] = mitjana(num,x)

mit_x = 0.*x;
if num == 1
    mit_x(1) = x(1);
    mit_x(2) = (x(1)+x(2))/2;
    mit_x(3) = (x(1)+x(2)+x(3))/3;
    mit_x(4) = (x(1)+x(2)+x(3)+x(4))/4;
    mit_x(5) = (x(1)+x(2)+x(3)+x(4)+x(5))/5;
    mit_x(6) = (x(1)+x(2)+x(3)+x(4)+x(5)+x(6))/6;
    for i = 7:length(x)
        mit_x(i) = (x(i)+x(i-1)+x(i-2)+x(i-3)+x(i-4)+x(i-5)+x(i-6))/7;
    end
end
if num == 2
    mit_x(1) = x(1);
    mit_x(2) = (x(1)+x(2))/2;
    mit_x(3) = (x(1)+x(2)+x(3))/3;
    for i = 4:length(x)-3
        mit_x(i) = (x(i)+x(i-1)+x(i-2)+x(i-3)+x(i+1)+x(i+2)+x(i+3))/7;
    end
    mit_x(end-2) = (x(end-1)+x(end-2)+x(end))/3;
    mit_x(end-1) = (x(end-1)+x(end))/2;
    mit_x(end) = x(end);
end
if num == 3
    mit_x(end) = x(end);
    mit_x(end-1) = (x(end-1)+x(end))/2;
    mit_x(end-2) = (x(end-1)+x(end-2)+x(end))/3;
    mit_x(end-3) = (x(end-1)+x(end-2)+x(end-3)+x(end))/4;
    mit_x(end-4) = (x(end-1)+x(end-2)+x(end-3)+x(end-4)+x(end))/5;
    mit_x(end-5) = (x(end-1)+x(end-2)+x(end-3)+x(end-4)+x(end-5)+x(end))/6;
    for i = 1:length(x)-6
        mit_x(i) = (x(i)+x(i+1)+x(i+2)+x(i+3)+x(i+4)+x(i+5)+x(i+6))/7;
    end
end

end
function [] = write_analysis_errors(filename,time,countries,ERRORS_7,ERRORS_14,ERRORS_21)

TI = cell(length(time),1);
for k = 1:length(time)
    TI{k} = datestr(time(k),'dd/mm/yyyy');
end
pagina = 'Hoja1';
xlswrite(filename,TI,pagina,'A2')
xlswrite(filename,countries,pagina,'B2')
xlswrite(filename,round(ERRORS_7,3),pagina,'C2')
xlswrite(filename,round(ERRORS_14,3),pagina,'Y2')
xlswrite(filename,round(ERRORS_21,3),pagina,'AU2')

xlswrite(filename,{'Date'},pagina,'A1')
xlswrite(filename,{'Country'},pagina,'B1')
xlswrite(filename,{'ER7 7'},pagina,'C1')
xlswrite(filename,{'ER7 8'},pagina,'D1')
xlswrite(filename,{'ER7 9'},pagina,'E1')
xlswrite(filename,{'ER7 10'},pagina,'F1')
xlswrite(filename,{'ER7 11'},pagina,'G1')
xlswrite(filename,{'ER7 12'},pagina,'H1')
xlswrite(filename,{'ER7 13'},pagina,'I1')
xlswrite(filename,{'ER7 14'},pagina,'J1')
xlswrite(filename,{'ER7 15'},pagina,'K1')
xlswrite(filename,{'ER7 16'},pagina,'L1')
xlswrite(filename,{'ER7 17'},pagina,'M1')
xlswrite(filename,{'ER7 18'},pagina,'N1')
xlswrite(filename,{'ER7 19'},pagina,'O1')
xlswrite(filename,{'ER7 20'},pagina,'P1')
xlswrite(filename,{'ER7 21'},pagina,'Q1')
xlswrite(filename,{'ER7 22'},pagina,'R1')
xlswrite(filename,{'ER7 23'},pagina,'S1')
xlswrite(filename,{'ER7 24'},pagina,'T1')
xlswrite(filename,{'ER7 25'},pagina,'U1')
xlswrite(filename,{'ER7 26'},pagina,'V1')
xlswrite(filename,{'ER7 27'},pagina,'W1')
xlswrite(filename,{'ER7 28'},pagina,'X1')
xlswrite(filename,{'ER14 7'},pagina,'Y1')
xlswrite(filename,{'ER14 8'},pagina,'Z1')
xlswrite(filename,{'ER14 9'},pagina,'AA1')
xlswrite(filename,{'ER14 10'},pagina,'AB1')
xlswrite(filename,{'ER14 11'},pagina,'AC1')
xlswrite(filename,{'ER14 12'},pagina,'AD1')
xlswrite(filename,{'ER14 13'},pagina,'AE1')
xlswrite(filename,{'ER14 14'},pagina,'AF1')
xlswrite(filename,{'ER14 15'},pagina,'AG1')
xlswrite(filename,{'ER14 16'},pagina,'AH1')
xlswrite(filename,{'ER14 17'},pagina,'AI1')
xlswrite(filename,{'ER14 18'},pagina,'AJ1')
xlswrite(filename,{'ER14 19'},pagina,'AK1')
xlswrite(filename,{'ER14 20'},pagina,'AL1')
xlswrite(filename,{'ER14 21'},pagina,'AM1')
xlswrite(filename,{'ER14 22'},pagina,'AN1')
xlswrite(filename,{'ER14 23'},pagina,'AO1')
xlswrite(filename,{'ER14 24'},pagina,'AP1')
xlswrite(filename,{'ER14 25'},pagina,'AQ1')
xlswrite(filename,{'ER14 26'},pagina,'AR1')
xlswrite(filename,{'ER14 27'},pagina,'AS1')
xlswrite(filename,{'ER14 28'},pagina,'AT1')
xlswrite(filename,{'ER21 7'},pagina,'AU1')
xlswrite(filename,{'ER21 8'},pagina,'AV1')
xlswrite(filename,{'ER21 9'},pagina,'AW1')
xlswrite(filename,{'ER21 10'},pagina,'AX1')
xlswrite(filename,{'ER21 11'},pagina,'AY1')
xlswrite(filename,{'ER21 12'},pagina,'AZ1')
xlswrite(filename,{'ER21 13'},pagina,'BA1')
xlswrite(filename,{'ER21 14'},pagina,'BB1')
xlswrite(filename,{'ER21 15'},pagina,'BC1')
xlswrite(filename,{'ER21 16'},pagina,'BD1')
xlswrite(filename,{'ER21 17'},pagina,'BE1')
xlswrite(filename,{'ER21 18'},pagina,'BF1')
xlswrite(filename,{'ER21 19'},pagina,'BG1')
xlswrite(filename,{'ER21 20'},pagina,'BH1')
xlswrite(filename,{'ER21 21'},pagina,'BI1')
xlswrite(filename,{'ER21 22'},pagina,'BJ1')
xlswrite(filename,{'ER21 23'},pagina,'BK1')
xlswrite(filename,{'ER21 24'},pagina,'BL1')
xlswrite(filename,{'ER21 25'},pagina,'BM1')
xlswrite(filename,{'ER21 26'},pagina,'BN1')
xlswrite(filename,{'ER21 27'},pagina,'BO1')
xlswrite(filename,{'ER21 28'},pagina,'BP1')

end
function [] = write_analysis_sum(filename,time,countries,suma_7,suma_14,suma_21,e14_14)

TI = cell(length(time),1);
for k = 1:length(time)
    TI{k} = datestr(time(k),'dd/mm/yyyy');
end
pagina = 'Hoja1';
xlswrite(filename,TI,pagina,'A2')
xlswrite(filename,countries,pagina,'B2')
xlswrite(filename,round(suma_7,3),pagina,'C2')
xlswrite(filename,round(suma_14,3),pagina,'Y2')
xlswrite(filename,round(suma_21,3),pagina,'AU2')
xlswrite(filename,round(e14_14,3),pagina,'BQ2')

xlswrite(filename,{'Date'},pagina,'A1')
xlswrite(filename,{'Country'},pagina,'B1')
xlswrite(filename,{'SM7 7'},pagina,'C1')
xlswrite(filename,{'SM7 8'},pagina,'D1')
xlswrite(filename,{'SM7 9'},pagina,'E1')
xlswrite(filename,{'SM7 10'},pagina,'F1')
xlswrite(filename,{'SM7 11'},pagina,'G1')
xlswrite(filename,{'SM7 12'},pagina,'H1')
xlswrite(filename,{'SM7 13'},pagina,'I1')
xlswrite(filename,{'SM7 14'},pagina,'J1')
xlswrite(filename,{'SM7 15'},pagina,'K1')
xlswrite(filename,{'SM7 16'},pagina,'L1')
xlswrite(filename,{'SM7 17'},pagina,'M1')
xlswrite(filename,{'SM7 18'},pagina,'N1')
xlswrite(filename,{'SM7 19'},pagina,'O1')
xlswrite(filename,{'SM7 20'},pagina,'P1')
xlswrite(filename,{'SM7 21'},pagina,'Q1')
xlswrite(filename,{'SM7 22'},pagina,'R1')
xlswrite(filename,{'SM7 23'},pagina,'S1')
xlswrite(filename,{'SM7 24'},pagina,'T1')
xlswrite(filename,{'SM7 25'},pagina,'U1')
xlswrite(filename,{'SM7 26'},pagina,'V1')
xlswrite(filename,{'SM7 27'},pagina,'W1')
xlswrite(filename,{'SM7 28'},pagina,'X1')
xlswrite(filename,{'SM14 7'},pagina,'Y1')
xlswrite(filename,{'SM14 8'},pagina,'Z1')
xlswrite(filename,{'SM14 9'},pagina,'AA1')
xlswrite(filename,{'SM14 10'},pagina,'AB1')
xlswrite(filename,{'SM14 11'},pagina,'AC1')
xlswrite(filename,{'SM14 12'},pagina,'AD1')
xlswrite(filename,{'SM14 13'},pagina,'AE1')
xlswrite(filename,{'SM14 14'},pagina,'AF1')
xlswrite(filename,{'SM14 15'},pagina,'AG1')
xlswrite(filename,{'SM14 16'},pagina,'AH1')
xlswrite(filename,{'SM14 17'},pagina,'AI1')
xlswrite(filename,{'SM14 18'},pagina,'AJ1')
xlswrite(filename,{'SM14 19'},pagina,'AK1')
xlswrite(filename,{'SM14 20'},pagina,'AL1')
xlswrite(filename,{'SM14 21'},pagina,'AM1')
xlswrite(filename,{'SM14 22'},pagina,'AN1')
xlswrite(filename,{'SM14 23'},pagina,'AO1')
xlswrite(filename,{'SM14 24'},pagina,'AP1')
xlswrite(filename,{'SM14 25'},pagina,'AQ1')
xlswrite(filename,{'SM14 26'},pagina,'AR1')
xlswrite(filename,{'SM14 27'},pagina,'AS1')
xlswrite(filename,{'SM14 28'},pagina,'AT1')
xlswrite(filename,{'SM21 7'},pagina,'AU1')
xlswrite(filename,{'SM21 8'},pagina,'AV1')
xlswrite(filename,{'SM21 9'},pagina,'AW1')
xlswrite(filename,{'SM21 10'},pagina,'AX1')
xlswrite(filename,{'SM21 11'},pagina,'AY1')
xlswrite(filename,{'SM21 12'},pagina,'AZ1')
xlswrite(filename,{'SM21 13'},pagina,'BA1')
xlswrite(filename,{'SM21 14'},pagina,'BB1')
xlswrite(filename,{'SM21 15'},pagina,'BC1')
xlswrite(filename,{'SM21 16'},pagina,'BD1')
xlswrite(filename,{'SM21 17'},pagina,'BE1')
xlswrite(filename,{'SM21 18'},pagina,'BF1')
xlswrite(filename,{'SM21 19'},pagina,'BG1')
xlswrite(filename,{'SM21 20'},pagina,'BH1')
xlswrite(filename,{'SM21 21'},pagina,'BI1')
xlswrite(filename,{'SM21 22'},pagina,'BJ1')
xlswrite(filename,{'SM21 23'},pagina,'BK1')
xlswrite(filename,{'SM21 24'},pagina,'BL1')
xlswrite(filename,{'SM21 25'},pagina,'BM1')
xlswrite(filename,{'SM21 26'},pagina,'BN1')
xlswrite(filename,{'SM21 27'},pagina,'BO1')
xlswrite(filename,{'SM21 28'},pagina,'BP1')
xlswrite(filename,{'E14 14'},pagina,'BQ1')

end
