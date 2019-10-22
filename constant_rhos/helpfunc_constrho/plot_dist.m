clear all
close all
clc

shape_alpha = [6];
scale_alpha = [0.008, 0.01, 0.02, 0.03, 0.05, 0.08];%, 0.1, 0.3];
alpha_max = 1;
xx = linspace(0, 1, 500);


k = 0;
figure
for i = 1:length(shape_alpha)
	for j = 1:length(scale_alpha)
		k = k + 1;
		pdf = gampdf(xx, shape_alpha(i), scale_alpha(j));
		leg{k} = ['\alpha = ', num2str(shape_alpha(i)), '; \beta = ', num2str(scale_alpha(j))];
		
		plot(xx, pdf, 'LineWidth', 2)
		hold on
	end
end
legend(leg)
xlabel('Random (1-porosity)')
ylabel('gamma PDF')


% Results
% ex: shape = 10, scale = 0.01 => 0.1 mode, 0.1 width
% mode = shape * scale (approx)
% Good set: shape = 6, scale = [0.008: 0.08]


%% Lognormal vs. normal
% 
% mean_alpha = 1;
% std_alpha = [0.01, 0.03, 0.05, 0.08];
% xx = linspace(0.5, 2, 500);
% 
% 
% k = 0;
% figure
% for i = 1:length(mean_alpha)
% 	for j = 1:length(std_alpha)
%         
% 		k = k + 1;
%         mu(k) = log(mean_alpha(i).^2./(sqrt(std_alpha(j) + mean_alpha(i).^2)));
%         sig(k) = sqrt(log(std_alpha(j)./mean_alpha(i).^2 + 1));
% 
% 
% 		pdfn = lognpdf(xx, mu(k), sig(k));
%         pdfg = normpdf(xx, mean_alpha(i), std_alpha(j));
% 		leg{2*k-1} = ['\alpha = ', num2str(mu(k)), '; \sigma = ', num2str(sig(k))];
% 		leg{2*k} = ['\alpha = ', num2str(mean_alpha(i)), '; \sigma = ', num2str(std_alpha(j))];
%         
% 		plot(xx, pdfn, '-', xx, pdfg, '--', 'LineWidth', 2)
% 		hold on
% 	end
% end
% legend(leg)
% title('comparison of normal pdf and lognormal pdf')
% xlabel('Random (porosity)')
% ylabel('Log noraml PDF')

%% Lognormal

mean_alpha = .2;
std_alpha = [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.09];
xx = linspace(0, 2, 500);


k = 0;
figure
for i = 1:length(mean_alpha)
	for j = 1:length(std_alpha)
        
		k = k + 1;
        mu(k) = log(mean_alpha(i).^2./(sqrt(std_alpha(j) + mean_alpha(i).^2)));
        sig(k) = sqrt(log(std_alpha(j)./mean_alpha(i).^2 + 1));

		pdfn = lognpdf(xx, mu(k), sig(k));
		leg{k} = ['\alpha = ', num2str(mu(k)), '; \sigma = ', num2str(sig(k))];
        
		plot(1 - xx, pdfn, '-', 'LineWidth', 2)
		hold on
	end
end

legend(leg)
title('Comparison of lognormal pdf')
xlabel('Random (porosity)')
ylabel('Log noraml PDF')

