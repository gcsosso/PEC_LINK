clear 
close all

%% Timeseries data
raw_data = readmatrix("Fig4b_Fluorescence measurements.xlsx");

Col0_f_data=raw_data(2:8,2:40);
Col0_t_data=[-3:3:15]';
Col0_I=mean(Col0_f_data,2);
Col0_I_var=var(Col0_f_data,0,2);
Col0_t=60*Col0_t_data;
Col0_t=Col0_t-Col0_t(1);

sfr8_f_data=raw_data(12:18,2:43);
sfr8_t_data=[-3:3:15]';
sfr8_I=mean(sfr8_f_data,2);
sfr8_I_var=var(sfr8_f_data,0,2);
sfr8_t=60*sfr8_t_data;
sfr8_t=sfr8_t-sfr8_t(1);

Col0_pore_area=1.0e-18*8.354;
sfr8_pore_area=1.0e-18*10.729;

Col0_pore_per_area=mean([930.6,996.3,1231.6])/(0.5e-6)^2;
sfr8_pore_per_area=mean([1016.8, 1195.0,1436.1])/(0.5e-6)^2;

Col0_R=Col0_pore_area*Col0_pore_per_area;
sfr8_R=sfr8_pore_area*sfr8_pore_per_area;



%% SV Data
conc_SV=[0.00237 0.50165 1.00083 1.99924]*1.0e-3;
sfr8_I_SV=1.0./[1.0283 4.5255 9.2897 20.0091];
Col0_I_SV=1.0./[0.4167 1.8452 4.5833 9.5238];
if(false)
figure
plot(conc_SV,1.0./sfr8_I_SV);
hold on
plot(conc_SV,1.0./Col0_I_SV);
end
%% Functions

I_func=@(q_m,K) 1./(1+K*q_m);
q_m_func=@(t,T,q_b)q_b*(1-exp(-t/T));


%% Fitting for K
q_b=(0.5e-6)*1.0e3;
ft=fittype('a*(1-exp(-x/T))');
 Col0_coeffs=fit(Col0_t,(1./Col0_I-1)/q_b,ft,'StartPoint',[1e4,1/Col0_t(end)]);
 sfr8_coeffs=fit(sfr8_t,(1./sfr8_I-1)/q_b,ft,'StartPoint',[1e4,1/Col0_t(end)]);
 a=coeffvalues(Col0_coeffs);
 Col0_K=a(2);
 Col0_T=a(1);
 a=coeffvalues(sfr8_coeffs);
 sfr8_K=a(2);
 sfr8_T=a(1);

 t_plot=linspace(0,3000,1000);
 fig=figure();
ax = gca();
set(fig,'Units','centimeters');
fig_aspect_ratio=9/16;
text_width=13.49414;
set(fig,'position',[2,2,text_width,fig_aspect_ratio*text_width]);
set(ax,'FontSize',9,'TickDir','out','TickLabelInterpreter','latex');
hold on

 cmap=brewermap(2,'Set1');
 errorbar(Col0_t/60-3,Col0_I,Col0_I_var,'o','DisplayName','Col0 data','LineWidth',2,'color',cmap(1,:),'MarkerSize',2)
 hold on
 errorbar(sfr8_t/60-3,sfr8_I,sfr8_I_var,'o','DisplayName','sfr8 data','LineWidth',2,'color',cmap(2,:),'MarkerSize',2)
 plot(t_plot/60-3,I_func(q_m_func(t_plot,Col0_T,q_b),Col0_K),'-','DisplayName','Col0 fit','LineWidth',2,'color',cmap(1,:))
  plot(t_plot/60-3,I_func(q_m_func(t_plot,sfr8_T,q_b),sfr8_K),'-','DisplayName','sfr8 fit','LineWidth',2,'color',cmap(2,:))
xlim([-3,30])
xlabel('Time (Minutes)')
ylabel('Relative Intensity')
  legend()
  exportgraphics(fig,'POSTREVIEW2_fitting.png','Resolution',300)
  exportgraphics(fig,'POSTREVIEW2_fitting.pdf')
    savefig('POSTREVIEW2_Fig3b.fig')

