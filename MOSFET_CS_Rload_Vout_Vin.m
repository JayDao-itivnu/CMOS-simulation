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
Vin1_n = subs(Vin1_syms,[Vdd Rd kn Vth],[10 1e3 1e-3 1.5]);% Symbolic substitution
Vin1 = double(Vin1_n(1)); % Take the positive result

%In saturation region: Vout = Vdd - Rd*kn*(Vin - Vth)^2

Vout_sat = Vdd - Rd*kn*(Vin - Vth)^2;
Vout_sat_n = subs(Vout_sat,[Vdd Rd kn Vth],[10 1e3 1e-3 1.5]);

%In triode region: Vout_tri = Vdd - Rd*kn*(2*(Vin-Vth)*Vout_tri-Vout_tri^2)

Vout_tri = solve(Vout_tri-(Vdd - Rd*kn*(2*(Vin-Vth)*Vout_tri-Vout_tri^2)),Vout_tri);
Vout_tri_n = subs(Vout_tri(2),[Vdd Rd kn Vth],[10 1e3 1e-3 1.5]);

% Input signal

Vin_0 = 0:0.1:10;
Vout_n = zeros(1,length(Vin_0));
gm_n = zeros(1,length(Vin_0));

for i=1:1:length(Vin_0)
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

hl1 = plot(Vin_0,Vout_n); % Input voltage
hl2 = plot([Vth_0 Vth_0],[0 11]);
hl3 = plot([Vin1 Vin1],[0 11]);

ax1 = gca;
set(ax1,'Xlim',[0 10]);
set(ax1,'Ylim',[0 11]);
set(ax1,'XColor','k','YColor','k');
set(get(ax1,'Title'),'String','Output - Input Voltage Characteristics','FontSize', 12);
set(get(ax1,'XLabel'),'String','Input Voltage - V','FontSize', 12);
set(get(ax1,'YLabel'),'String','Output Voltage - V','FontSize', 12);
set(ax1,'FontSize', 12);
set(ax1,'Box','On');

set(hl1,'LineWidth',2.5);
set(hl1,'LineStyle','-');
set(hl1,'Color','b');

set(hl2,'LineWidth',2);
set(hl2,'LineStyle','--');
set(hl2,'Color','k');

set(hl3,'LineWidth',2);
set(hl3,'LineStyle','--');
set(hl3,'Color','k');

text(Vth_0,0.5,'Vth');
text(Vin1,0.5,'Vin1');

% Id - Vin Characteristics

Id_n = (Vdd_0-Vout_n)/Rd_0;
%gm_n = sqrt(kn_0*Id_n);

figure(2), grid on, hold on,

hl1 = plot(Vin_0,Id_n); % Example 3.1(a) (page 49)
hl2 = plot([Vth_0 Vth_0],[0 0.01]);
hl3 = plot([Vin1 Vin1],[0 0.01]);

ax1 = gca;
set(ax1,'Xlim',[0 10]);
%set(ax1,'Ylim',[0 11]);
set(ax1,'XColor','k','YColor','k');
set(get(ax1,'Title'),'String','Drain current - Gate-Source Voltage Characteristics','FontSize', 12);
set(get(ax1,'XLabel'),'String','Input Voltage - V','FontSize', 12);
set(get(ax1,'YLabel'),'String','Drain Current - A','FontSize', 12);
set(ax1,'FontSize', 12);
set(ax1,'Box','On');

set(hl1,'LineWidth',2.5);
set(hl1,'LineStyle','-');
set(hl1,'Color','b');

set(hl2,'LineWidth',2);
set(hl2,'LineStyle','--');
set(hl2,'Color','k');

set(hl3,'LineWidth',2);
set(hl3,'LineStyle','--');
set(hl3,'Color','k');

text(Vth_0,0.001,'Vth');
text(Vin1,0.006,'Vin1');

figure(3), grid on, hold on,

hl1 = plot(Vin_0,gm_n); % Example 3.1(b) (page 49)

ax1 = gca;
set(ax1,'Xlim',[0 10]);
%set(ax1,'Ylim',[0 11]);
set(ax1,'XColor','k','YColor','k');
set(get(ax1,'Title'),'String','Transconductance - Input Voltage Characteristics','FontSize', 12);
set(get(ax1,'XLabel'),'String','Input Voltage - V','FontSize', 12);
set(get(ax1,'YLabel'),'String','Gm - S','FontSize', 12);
set(ax1,'FontSize', 12);
set(ax1,'Box','On');

set(hl1,'LineWidth',2.5);
set(hl1,'LineStyle','-');
set(hl1,'Color','b');

text(0.2,0.2e-3,'Turn-off');
text(2.5,1.2e-3,'Saturation');
text(6.5,2.5e-3,'Triode');

