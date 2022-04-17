%VNU.UET.FET.MEMS
%MOSFET
%Common-Source Stage Amplifier
%with Resistance load

clear all, close all,
syms Vin Vin1_syms Vout_sat Vout_sat_n Vout_tri Vout_tri_n kn Vth Rd Vdd

%MOS parameters:% kn = 1e-3; % kn = 1/2*umn*Cox*W/L, Vth = 1.5; % Threshold volgate
Vth_0 = 1.5;
kn_0 = 1e-3;
%Circuit parameters % Rd = 1e3; % Vdd = 10;
Vdd_0 = 10;
Rd_0 = 1e3;

% Solving for Vin1
Vin1_syms = solve((Vin1_syms-Vth)-(Vdd-Rd*kn*(Vin1_syms-Vth)^2),Vin1_syms);
Vin1_n = subs(Vin1_syms,[Vdd Rd kn Vth],[10 Rd_0 1e-3 1.5]);% Symbolic substitution
Vin1 = double(Vin1_n(1)); % Take the positive result

%In saturation region: Vout = Vdd - Rd*kn*(Vin - Vth)^2
Vout_sat = Vdd - Rd*kn*(Vin - Vth)^2;
Vout_sat_n = subs(Vout_sat,[Vdd Rd kn Vth],[10 Rd_0 1e-3 1.5]);

%In triode region: Vout_tri = Vdd - Rd*kn*(2*(Vin-Vth)*Vout_tri-Vout_tri^2)
Vout_tri = solve(Vout_tri-(Vdd - Rd*kn*(2*(Vin-Vth)*Vout_tri-Vout_tri^2)),Vout_tri);
Vout_tri_n = subs(Vout_tri(2),[Vdd Rd kn Vth],[10 Rd_0 1e-3 1.5]);

% Input signal
t = 0:1e-5:2e-3;
Vin_dc = 3.0; % Sinh vien thu voi cac dien ap bias khac nhau
Vin_am = 0.5; % Sinh vien thu voi cac bien do tin hieu loi vao khac nhau
Vin_f = 1e3;
Vin_0 = Vin_dc + Vin_am*cos(2*pi*Vin_f*t+3*pi/2);

Vout_n = zeros(1,length(t));
gm_n = zeros(1,length(Vin_0));

for i=1:1:length(t)
    if Vin_0(i) <= Vth_0 % Threshold voltage = 1.5 V
        Vout_n(i) = Vdd_0; % Turnoff 
        gm_n(i) = 0;
    elseif Vin_0(i) <= Vin1 % saturation region
        Vout_n(i) = double(subs(Vout_sat_n,Vin,Vin_0(i)));
        gm_n(i) = 2*kn_0*(Vin_0(i)-Vth_0);
    elseif Vin_0(i) > Vin1 % Triode region
        Vout_n(i) = double(subs(Vout_tri_n,Vin,Vin_0(i)));
        gm_n(i) = 2*kn_0*Vout_n(i); % Vds = Vout - in this case 
    end
end

figure(1), grid on, hold on,

hl1 = plot(t,Vin_0); % Input voltage
hl2 = plot(t,ones(1,length(t))*Vth_0); % Threshold voltage
hl3 = plot(t,ones(1,length(t))*Vin1); % Saturation - Triode transision voltage
hl4 = plot(t,Vout_n); % Output voltage

ax1 = gca;
set(ax1,'Xlim',[0 2e-3]);
set(ax1,'Ylim',[0 10]);
set(ax1,'XColor','k','YColor','k');
set(get(ax1,'XLabel'),'String','Time - s','FontSize', 14);
set(get(ax1,'YLabel'),'String','Voltage - V','FontSize', 14);
set(ax1,'FontSize', 12);
set(ax1,'Box','On');

set(hl1,'LineWidth',2);
set(hl1,'LineStyle','-');
set(hl1,'Color','r');

set(hl2,'LineWidth',2);
set(hl2,'LineStyle','--');
set(hl2,'Color','k');

set(hl3,'LineWidth',2);
set(hl3,'LineStyle','--');
set(hl3,'Color','m');

set(hl4,'LineWidth',2);
set(hl4,'LineStyle','-');
set(hl4,'Color','b');

legend('Vin',...
       'Vth',...
       'Vin1',...
       'Vout',...
       'Location','southeast')
   
figure(2), grid on, hold on,
subplot(3,1,1), grid on, hold on,
hl1 = plot(t,Vin_0); % Input voltage
hl2 = plot(t,ones(1,length(t))*Vth_0); % Threshold voltage
hl3 = plot(t,ones(1,length(t))*Vin1); % Saturation - Triode transision voltage
ylabel('Vin','FontSize', 12);

subplot(3,1,2), grid on, hold on,
hl4 = plot(t,Vout_n); % Output voltage
ylabel('Vout','FontSize', 12);

subplot(3,1,3), grid on, hold on,
hl5 = plot(t,gm_n);

ax1 = gca;
set(ax1,'XColor','k','YColor','k');
set(get(ax1,'XLabel'),'String','Time - s','FontSize', 12);
set(get(ax1,'YLabel'),'String','Gm','FontSize', 12);
set(ax1,'FontSize', 12);
set(ax1,'Box','On');

set(hl1,'LineWidth',2);
set(hl1,'LineStyle','-');
set(hl1,'Color','r');

set(hl2,'LineWidth',2);
set(hl2,'LineStyle','--');
set(hl2,'Color','k');

set(hl3,'LineWidth',2);
set(hl3,'LineStyle','--');
set(hl3,'Color','m');

set(hl4,'LineWidth',2);
set(hl4,'LineStyle','-');
set(hl4,'Color','b');

set(hl5,'LineWidth',2);
set(hl5,'LineStyle','-');
set(hl5,'Color','b');
