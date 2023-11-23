clear;
load altitude.mat

%% FP - AOA
figure(1)
plot(T(1:935), alt(1:935), "LineWidth", 3)
hold on
xlabel("Time (s)")
ylabel("Altitude (km)")

rp = plot(T(46), alt(46), '^', "Color", "#FF8800", "MarkerSize", 12);
meco = plot(T(251), alt(251), '^', "Color", "#FF8800", "MarkerSize", 12);
rp_fin = plot(T(239), alt(239), '^', "Color", "#FF8800", "MarkerSize", 12);
coasting = plot(T(283), alt(283), '^', "Color", "#FF8800", "MarkerSize", 12);
circ_burn = plot(T(922), alt(922), '^', "Color", "#FF8800", "MarkerSize", 12);
burn_over = plot(T(933), alt(933), '^', "Color", "#FF8800", "MarkerSize", 12);


text(T(46) - 15, alt(46) + 50, "Roll program start "+ newline + " T+22.5 s")
text(T(239) - 25 , alt(239) - 45, "Roll program end" + newline + "T+215 s")
text(T(251) + 30, alt(251) - 10, "MECO & Hot staging" + newline + " T+227 s")
text(T(283) + 50, alt(283) + 20, "Coasting" + newline + "T+243 s")
text(T(922) - 150, alt(922)+ 10, "Circularization burn"+ newline + " T+633 s")
text(T(933) + 20, alt(933)- 20, "SECO"+ newline + " T+647 s")


ax = gca;
ax.FontSize = 12;
set(gcf, 'Position',  [600, 400, 1000, 600])
grid on
