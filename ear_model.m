clearvars
close all

function [num, den] = polar_to_cartesian(zero, pole)
    z = numel(zero);
    p = numel(pole);
    num = 1;
    den = 1;
    for j = 1:z
        num = conv(num, [-zero(j), 1]);
    end
    for j = 1:p
        den = conv(den, [-pole(j), 1]);
    end
end

function h_z_arr = transfer_function(num, den, arr)
    n = numel(num);
    d = numel(den);
    l = numel(arr);
    arr_pow = ones(1, l);
    num_sum = zeros(1, l);
    den_sum = zeros(1, l);
    for j = 1:max(n, d)
        if j <= n
            num_iter = num(j) .* arr_pow;
            num_sum = num_sum + num_iter; 
        end
        if j <= d
            den_iter = den(j) .* arr_pow;
            den_sum = den_sum + den_iter;
        end
        arr_pow = arr_pow .* arr;
    end
    h_z_arr = num_sum ./ den_sum;
end

%% Outer Ear---------------------------------------------------------------

z_outer = [0.8*exp(2i*pi*2/200), 0.8*exp(-2i*pi*2/200), 0.8*exp(2i*pi*95/200), 0.8*exp(-2i*pi*95/200)];
p_outer = [0, 0.9*exp(2i*pi*3/20), 0.9*exp(-2i*pi*3/20), 0.6*exp(2i*pi*5/20), 0.6*exp(-2i*pi*5/20)];

[n_outer, d_outer] = polar_to_cartesian(z_outer, p_outer);
n_outer = n_outer.*10;

f = (100:100:10000);
exps = exp((2i*pi/20000) .* f);

freq_resp_outer = transfer_function(n_outer, d_outer, exps);
mag_resp_outer = mag2db(abs(10.*freq_resp_outer));

figure
plot(f, mag_resp_outer);
set(gca, XScale='log');
title('Approx. magnitude response of outer ear');

%% Middle Ear--------------------------------------------------------------

z_middle = [
    0.99*exp(2i*pi*20/20000), 0.99*exp(-2i*pi*20/20000), -0.1, 0.7...
    0.5*exp(2i*pi*2000/20000), 0.5*exp(-2i*pi*2000/20000),...
    ];
p_middle = [
    0, 0.8*exp(2i*pi*200/20000), 0.8*exp(-2i*pi*200/20000),...
    0.9*exp(2i*pi*800/20000), 0.9*exp(-2i*pi*800/20000), ...
    0.7*exp(2i*pi*3000/20000), 0.7*exp(-2i*pi*3000/20000),...
    ];

[n_middle, d_middle] = polar_to_cartesian(z_middle, p_middle);
n_middle = n_middle.*(10^1.25);

f = (100:100:10000);
exps = exp((2i*pi/20000) .* f);

freq_resp_middle = transfer_function(n_middle, d_middle, exps);
mag_resp_middle = mag2db(abs(freq_resp_middle));

figure
plot(f, mag_resp_middle);
set(gca, XScale='log');
title('Approx. magnitude response of middle ear');

%% Combined response-------------------------------------------------------

n_combined = conv(n_outer, n_middle);
d_combined = conv(d_outer, d_middle);

f = (100:100:10000);
exps = exp((2i*pi/20000) .* f);

freq_resp_combined = transfer_function(n_combined, d_combined, exps);
mag_resp_combined = mag2db(abs(freq_resp_combined));

figure
plot(f, mag_resp_combined);
set(gca, XScale='log');
title('Approx. magnitude response of outer+middle ear');

%% Inner ear (Cochlea): Constant setup-------------------------------------

T = 1/48000;
filter_count = 128;
filter_num = (1:filter_count)' .* 3.5/filter_count;

f_p = (20000) .* (10 .^ (-2/3 .* filter_num));
f_z = f_p .* 1.05;
Q_p = linspace(10, 5.5, filter_count)';
Q_z = linspace(22, 12, filter_count)';
K = f_p ./ f_z;

% p_1 = pi .* (f_z ./ Q_z);
% q_1 = p_1 .* (sqrt(4.*(Q_z.^2) - 1));
% p_2 = pi .* (f_p ./ Q_p);
% q_2 = p_2 .* (sqrt(4.*(Q_p.^2) - 1));

theta_z = cos((2*pi*T).*f_z);
r_z = exp(-(pi*T).*(f_z./Q_z));
a_1 = 2 .* r_z .* theta_z;
a_2 = r_z.^2;

theta_p = cos((2*pi*T).*f_p);
r_p = exp(-(pi*T).*(f_p./Q_p));
b_1 = 2 .* r_p .* theta_p;
b_2 = r_p.^2;

theta_lp = cos((2*pi*T*sqrt(2)).*f_z);
a_0 = 2 - theta_lp - sqrt((2 - theta_lp).^2 - 1);
% a_1 = 2 .* exp(-T .* p_1) .* cos(T .* q_1);
% a_2 = exp(-2*T .* p_1);
% b_1 = 2 .* exp(-T .* p_2) .* cos(T .* q_2);
% b_2 = exp(-2*T .* p_2);
gain = K .* (1 - a_0) .* (1 - b_1 + b_2);

pressure(filter_count) = struct('num', [], 'den', []);
displacement(filter_count) = struct('num', [], 'den', []);

for j = 1:filter_count
    den_iter = conv([1, -a_0(j)], [1, -b_1(j), b_2(j)]);

    pressure(j).num = gain(j) .* [1, -a_1(j), a_2(j)];
    pressure(j).den = den_iter .* (1-a_1(j)+a_2(j));

    displacement(j).num = gain(j);
    displacement(j).den = den_iter;
end

%% Inner ear (Cochlea): Trial and error------------------------------------

function [y_d, y_p] = displacement_output(filter_no, pres, disp, x)
    y_iter = x;
    if (filter_no > 1)
        for j = 1:filter_no-1
            y_iter = filter(pres(j).num, pres(j).den, y_iter);
        end
    end
    y_d = filter(disp(filter_no).num, disp(filter_no).den, y_iter);
    y_p = filter(pres(filter_no).num, pres(filter_no).den, y_iter);
end

imp_len = 8000;
impulse = [1, zeros(1, imp_len-1)];
plot_filters = [30, 60, 90];

figure
for j = 1:numel(plot_filters)
    h_z_1 = displacement_output(plot_filters(j), pressure, displacement, impulse);
    h_z_2 = displacement_output(plot_filters(j)-1, pressure, displacement, impulse);
    subplot(numel(plot_filters), 1, j)
    plot(h_z_1 - h_z_2);
    title(['Impulse response; f_p = ' num2str(f_p(plot_filters(j)))], FontName= 'Palatino Linotype')
end

sinusoid = sin((2*pi*3200) .* linspace(0, 1, 1/T))...
            + 1*sin((2*pi*1000) .* linspace(0, 1, 1/T))...
            + 1*sin((2*pi*300) .* linspace(0, 1, 1/T));

out_samples = (2000:50:8000);
freqHz = (0:filter_count-3)/(2*T*filter_count);
sin_d = zeros(filter_count, length(sinusoid));
sin_p = zeros(filter_count, length(sinusoid));
for j = 1:filter_count
    [sin_d(j, :), sin_p(j, :)] = displacement_output(j, pressure, displacement, sinusoid);
end

sin_d = diff(sin_d, 1, 1);
sin_d_prime = sin_d(2:end, :)-sin_d(1:end-1, :);

figure
subplot(2, 1, 1)
plot(sin_p(:, out_samples), Color= '#000000');
title('Pressure accross ear canal', FontName='Palatino Linotype');
subplot(2, 1, 2)
plot(sin_d_prime(:, out_samples), Color= '#000000');
title('Differentiated displacement accross ear canal', FontName='Palatino Linotype');

f = [(100:10:1990), (2000:100:24000)];
exps = exp((2i*pi*T) .* f);
freq_resp_inner = ones(1, numel(f));
filter_no = 48;
for j = 1:filter_no-1
    freq_resp_inner = freq_resp_inner .* transfer_function(pressure(j).num, pressure(j).den, exps);
end
freq_resp_inner = freq_resp_inner .* transfer_function(displacement(filter_no).num, displacement(filter_no).den, exps);
mag_resp_inner = mag2db(abs(freq_resp_inner));

figure
plot(f, mag_resp_inner);
set(gca, XScale='log');
title(['Magnitude response of filter output w/  f_p = ' num2str(f_p(filter_no))], FontName='Palatino Linotype');

filter_no = 81;
freq_resp_disp = transfer_function(displacement(filter_no).num, displacement(filter_no).den, exps);
freq_resp_pres = transfer_function(pressure(filter_no).num, pressure(filter_no).den, exps);
mag_resp_disp = mag2db(abs(freq_resp_disp));
mag_resp_pres = mag2db(abs(freq_resp_pres));

figure
plot(f, mag_resp_disp);
hold on
plot(f, mag_resp_pres);
hold off
set(gca, XScale='log');
title(['Filter w/ f_p = ' num2str(f_p(filter_no)) ';   f_z = ' num2str(f_z(filter_no))], FontName= 'Palatino Linotype');

c_0 = exp(-2*pi*30*T);
sin_d_prime_relu = max(0, sin_d_prime);
sin_d_energy = zeros(size(sin_d_prime, 1), numel(sinusoid));

sin_d_energy(:, 1) = (1-c_0) .* sin_d_prime_relu(:, 1);
for j = 2:size(sin_d_energy, 2)
    sin_d_energy(:, j) = (1-c_0)*sin_d_prime_relu(:, j) + c_0*sin_d_energy(:, j-1);
end

figure
subplot(2, 1, 1)
plot(sin_d_prime(:, 4000));
subplot(2, 1, 2)
plot(sin_d_energy(:, 4000));

%% Goofing around