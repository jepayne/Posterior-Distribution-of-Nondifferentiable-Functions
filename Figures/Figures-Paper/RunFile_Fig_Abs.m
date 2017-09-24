% ----------------------------------------------------------------------------
% Absolute Value Example
% ----------------------------------------------------------------------------
% Description:  This file comptutes the plots for Figures 1 and 2 in the
%   paper Kitagawa, Montiel-Olea, Payne. "Posterior Distribution of
%   Non-Differential Functions". (2017).
% ----------------------------------------------------------------------------
% Authors:      Toru Kitagawa (t.kitagawa@ucl.ac.uk),
%               Jose-Luis Montiel-Olea (montiel.olea@nyu.edu), and
%               Jonathan Payne (jep459@nyu.edu)
% ----------------------------------------------------------------------------

%% ------------------------------------------------------------
% 0. Set up Parameters
% ------------------------------------------------------------
clear; clc;
rng('default') % Set random seed to default

n = 100; % Sample size

% Set size parameters for quantile:
m = 500; % Number of samples generated (either as bootstrap draws or posterior draws) [Set to 1000 in paper]
nZ = 500; % Number of draws to estimate parametric bootstrap [Set to 2000 in paper]
I = 500; % Number of draws from MCMC when using random MLE [Set to 2000 in paper]
IG = 10000; % Number of draws from MCMC when using MLE grid [Set to 10000 in paper]
IB = 500; %Number of bootstrap repetitions [Set to 2000 in paper]
sthetaML=-3:.5:3; % thetaML: Range of scaled MLE values (n^(1/2)thetaML) at which posterior is evaluated
lenML = size(sthetaML,2);

% Set size parameters for coverage calculations:
J = 101; % Number of grid points for theta in coverage plots
width = 2;
theta = linspace(-width,width,J)'; % Grid for the parameter space for coverage plot
mcov = 200; % Number of data sets for coverage [Set to 750 in paper]
Icov = 200; % Number of draws from MCMC for coverage [Set to 750 in paper]

ifig = 1; % Start figure count

% True model from which the data is generated
%   X_i ~ N(\theta0,\sig0)
theta0 = 0; % True value of theta from which the data is generated
sig0 = 1; % True value of theta from which the data is generated

% Transformation
g=@(x)(abs(x));
gprime=@(x,the)((the ~= 0)*sign(the)*x + (the == 0)*abs(x));

% Generate data
Z=randn(nZ,1); % Bootstrap draws and posterior draws (for analytical posterior)
dataMC = randn(m,n)+theta0; % Posterior draws (for MCMC models)

% Generate the quantiles of the parametric bootstrap for fixed sthetML
approx=abs(bsxfun(@plus,Z/(n^.5),(n^(-.5))*sthetaML));
Approx_quant_2S=quantile(approx,[.975,.025],1); % Two sided approximate CS
Approx_quant_1S=quantile(approx,[.975,0],1); % One sided approximate CS

% Indicator variables (to indicate which sections of the code are run)
quantnorm = 1; % If 1, then calculate the quantile plots for Normal priors
covnorm = 1; % If 1, then calculate the coverage plots for Normal priors
quantgamma_a = 1; % If 1, then calculate the quantile plots for Gamma priors using random MLE
quantgamma_b = 1; % If 1, then calculate the quantile plots for Gamma priors using MLE grid
covgamma = 1; % If 1, then calculate the coverage plots for Gamma priors
quantbeta_a = 1; % If 1, then calculate the quantile plots for Beta priors using random ML
quantbeta_b = 1; % If 1, then calculate the quantile plots for Beta priors using MLE grid
covbeta = 1; % If 1, then calculate the coverage plots for Beta priors

%% ------------------------------------------------------------
% 1. Figure: Normal priors with lambda2 = 5, 10
% ------------------------------------------------------------
disp('Normal Prior: Theta ~ N(0,1/lambda^2)')

for lambda2 = [5,10]
    % 1.1.) Calculate quantiles for a set vector of data realisation:
    % ------------------------------------------------------------
    if quantnorm == 1
        % 1.1.1.) Generate posterior draws
        % Model X_i ~ N(\theta,1) with prior N(0,(1/lambda2))
        post=abs(bsxfun(@plus,Z/((n+lambda2)^.5),((n^.5)/(n+lambda2))*sthetaML));

        % 1.1.2.) Calculate one sided and two sided confidence sets
        BCset1_quant=quantile(post,[0.95,0],1); % One sided confidence interval
        BCset2_quant=quantile(post,[.975,.025],1); % Two sided confidence interval

        % 1.1.3.) Plot two sided confidence set
        figure(ifig);
        step=sthetaML(1,end)-sthetaML(1,end-1);
        for ix=1:size(sthetaML,2)
            plot([sthetaML(1,ix)  sthetaML(1,ix)],[BCset2_quant(2,ix) BCset2_quant(1,ix)],'-.red'); hold on
            plot([sthetaML(1,ix)  sthetaML(1,ix)],[Approx_quant_2S(2,ix) Approx_quant_2S(1,ix)],':blue'); hold on
            plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[BCset2_quant(1,ix) BCset2_quant(1,ix)],'red'); hold on
            plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[BCset2_quant(2,ix) BCset2_quant(2,ix)],'red'); hold on
            plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(1,ix) Approx_quant_2S(1,ix)],'blue'); hold on
            plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(2,ix) Approx_quant_2S(2,ix)],'blue'); hold on

        end
        legend({'95% Credible Set based on the posterior quantiles','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11);
        legend('boxoff')
        axis([-3.5 3.5 0 1])
        xlabel('$Z_n = \sqrt{n} \hat{\theta}_n$','Interpreter','latex');
        title(['Credible Sets and Bootstrap CIs for $|\theta|$ with Prior $\theta \sim N(0,1/\lambda^2)$, $\lambda^2=$' num2str(lambda2)],'Interpreter','latex');
        hold off
        saveTightFigure(figure(ifig), strcat('abs1_','n=',num2str(n),'lambda2=',num2str(lambda2),'_95CS_2Side.pdf'));
        ifig = ifig + 1;

    end


    % 1.2.) Calculate coverage for a set of thetas:
    % ------------------------------------------------------------
    if covnorm == 1
        % 1.2.1.) Calculate the set of data realisations
        thetaML = bsxfun(@plus,(1/(n^.5))*randn(1,I),theta); % MLE from data realisations (This is J x I)

        % 1.2.2.) Set up containers for quantiles
        cposthigh = zeros(J, I);
        cpostlow = zeros(J, I);
        cboothigh = zeros(J, I);
        cbootlow = zeros(J, I);

        % 1.2.3.) Calculate quantiles
        for i = 1:J
            i
            % Calculate quantiles for the posterior
            post = g(bsxfun(@plus,Z'/((n+lambda2)^.5), (n/(n+lambda2))*reshape(thetaML(i,:), [I,1])));
            cposthigh(i, :) = reshape(quantile(post, .975,2),[1,I]);
            cpostlow(i, :) = reshape(quantile(post, .0275,2),[1,I]);
            % Calculate quantiles for the bootstrap
            boot = g(bsxfun(@plus,Z'/(n^.5), reshape(thetaML(i,:), [I,1])));
            cboothigh(i, :) = reshape(quantile(boot, .975,2),[1,I]);
            cbootlow(i, :) = reshape(quantile(boot, .0275,2),[1,I]);
        end
        coveragepost = mean(bsxfun(@le,cpostlow,g(theta)).*...
            bsxfun(@le,g(theta),cposthigh),2);
        coverageboot = mean(bsxfun(@le,cbootlow,g(theta)).*...
            bsxfun(@le,g(theta),cboothigh),2);

        % 1.2.4.) Report the coverage of this procedure over theta at deltstar
        figure(ifig);
        plot(theta,coveragepost,'-.red'); hold on
        plot(theta,coverageboot,'-.blue'); hold on
        axis([-2 2 0 1.2])
        legend({'95% Credible Set based on the posterior quantiles','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11);
        legend('boxoff')
        xlabel('\theta')
        ylabel('Coverage')
        title(['Coverage Probability of 95\% CS and Bootstrap CI for $|\theta|$ with Prior $\theta \sim N(0,1/\lambda^2)$, $\lambda^2=$' num2str(lambda2)],'Interpreter','latex');
        hold off
        saveTightFigure(figure(ifig), strcat('abs2_','n=',num2str(n),'lambda2=',num2str(lambda2),'_95CS_Coverage.pdf'));
        ifig = ifig + 1;
    end

end


%% ------------------------------------------------------------
% 2. Figures: Gamma Prior
% ------------------------------------------------------------
disp('Gamma Prior: Theta ~ Z - (width + 1), Z ~ Gamma(alphap,betap)')

% 2.1. Set up functions:
% ------------------------------------------------------------
% tic;
% Prior Theta = Z - (width + 1), Z ~ Gamma(alphap,betap)
% Note that the algorithm can be adapted to a different prior by just
% changing definiton of the prior below and the function logpdfaux
alphap = 2;
betap = 2;
logpdfaux = @(x,data, sig0, alphap, betap) (-normlike([x,sig0],data) + log(gampdf(x+width+1,alphap, betap)));
logproppdf = @(x,y) -.5*((x-y).^2);
proprnd= @(x) x+randn;


% 2.2a. Calculate quantiles for a set vector of data realisation (using random samples):
% ------------------------------------------------------------
if quantgamma_a == 1
    % 2.2a.1.) Set up the containers
    thetadata=zeros(m,1);
    qgthetahigh_2S=zeros(m,1); % Container for upper bound on 2 sided CS
    qgthetalow_2S=zeros(m,1); % Container for lower bound on 2 sided CS

    for i=1:m
        i
        % Generate I posterior draws from posterior distribution
        logpdf=@(x) logpdfaux(x,dataMC(i,:),sig0, alphap, betap);
        postMCMC = mhsample (0, I, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',100);
        % Calculate the MLE for the data
        thetadata(i,:)=(n^.5)*mean(dataMC(i,:));
        % Calculate the quantiles
        qgthetahigh_2S(i,:)=quantile(abs(postMCMC),.975); % upper bound on 2 sided CS
        qgthetalow_2S(i,:)=quantile(abs(postMCMC),.025); % lower bound on 2 sided CS
        clear logpdf postMCMC
    end

    % 2.2a.2.) Bootstrap and Asymptotic Approximation (Two Sided)
    figure(ifig)
    step=sthetaML(1,end)-sthetaML(1,end-1);
    h1=scatter(thetadata,qgthetahigh_2S,1.5,'red'); hold on
    scatter(thetadata,qgthetalow_2S,1.5,'red'); hold on
    for ix=1:size(sthetaML,2)
    	h2=plot([sthetaML(1,ix)  sthetaML(1,ix)],[Approx_quant_2S(2,ix) Approx_quant_2S(1,ix)],':blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(2,ix) Approx_quant_2S(2,ix)],'blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(1,ix) Approx_quant_2S(1,ix)],'blue'); hold on
    end
    axis([-3.5 3.5 0 1])
    legend([h1,h2],{'End Points of 95% CS from posterior quantiles (MCMC)','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11)
    legend('boxoff')
    box on
    hold off
    xlabel('$Z_n = \sqrt{n} \hat{\theta}_n$','Interpreter','latex');
    title(['Credible Sets and Bootstrap CIs for $|\theta|$ with Prior $\theta \sim \gamma(2,2)-3$'],'Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('abs1_','n=',num2str(n),'Posterior95csGamma_2Side_Rand.pdf'));
    ifig = ifig + 1;
end

% 2.2b. Calculate quantiles for a set vector of data realisation (using grid of MLE):
% ------------------------------------------------------------
normlike_explicit = @(x, thetaML, sig0) ( normpdf(x, thetaML, sig0/sqrt(n)) );
logpdfaux_explicit = @(x,thetaML, sig0, alphap, betap) (log(normlike_explicit(x, thetaML, sig0)) + log(gampdf(x+width+1,alphap, betap)));

if quantgamma_b == 1
    % 2.2b.1.) Set up the containers
    qgthetahigh_2S=zeros(lenML,1); % Container for upper bound on 2 sided CS
    qgthetalow_2S=zeros(lenML,1); % Container for lower bound on 2 sided CS

    for i=1:lenML
        i
        % Generate posterior draws
        logpdf=@(x) logpdfaux_explicit(x,sthetaML(i)/sqrt(n),sig0, alphap, betap);
        postMCMC = mhsample (0, IG, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',100);
        
        % Calculate the quantiles
        qgthetahigh_2S(i,:)=quantile(abs(postMCMC),.975); % upper bound on 2 sided CS
        qgthetalow_2S(i,:)=quantile(abs(postMCMC),.025); % lower bound on 2 sided CS
        clear logpdf postMCMC
    end

    % 2.2b.2.) Bootstrap and Asymptotic Approximation (Two Sided)
    figure(ifig)
    step=sthetaML(1,end)-sthetaML(1,end-1);
    for ix=1:size(sthetaML,2)
        h1=plot([sthetaML(1,ix)  sthetaML(1,ix)],[qgthetalow_2S(ix) qgthetahigh_2S(ix)],':red'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[qgthetahigh_2S(ix) qgthetahigh_2S(ix)],'red'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[qgthetalow_2S(ix) qgthetalow_2S(ix)],'red'); hold on
        
    	h2=plot([sthetaML(1,ix)  sthetaML(1,ix)],[Approx_quant_2S(2,ix) Approx_quant_2S(1,ix)],':blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(2,ix) Approx_quant_2S(2,ix)],'blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(1,ix) Approx_quant_2S(1,ix)],'blue'); hold on
    end
    axis([-3.5 3.5 0 1])
    legend([h1,h2],{'End Points of 95% CS from posterior quantiles (MCMC)','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11)
    legend('boxoff')
    box on
    hold off
    xlabel('$Z_n = \sqrt{n} \hat{\theta}_n$','Interpreter','latex');
    title(['Credible Sets and Bootstrap CIs for $|\theta|$ with Prior $\theta \sim \gamma(2,2)-3$'],'Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('abs1_','n=',num2str(n),'Posterior95csGamma_2Side.pdf'));
    ifig = ifig + 1;
end


% 2.3. Calculate coverage for a set vector of data realisation:
% ------------------------------------------------------------
if covgamma == 1
    % 2.3.1.) Calculate data realisations
    database = randn(mcov,n);
    thetaML = bsxfun(@plus,(1/(n^.5))*randn(1,mcov),theta); % MLE from data realisations (This is J x I)

    % % 2.3.2.) Set up containers for quantiles
    cposthigh = zeros(J, mcov);
    cpostlow = zeros(J, mcov);
    cboothigh = zeros(J, mcov);
    cbootlow = zeros(J, mcov);

    % 2.3.3.) Calculate quantiles
    for j = 1:J
        j
        dataMCtheta = database + theta(j);
        for i=1:mcov
            % Generate I posterior draws from posterior distribution
            logpdf = @(x) logpdfaux(x,dataMCtheta(i,:),sig0, alphap, betap);
            postMCMC = mhsample(0, Icov, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',100);
            % Calculate the quantiles
            cposthigh(j,i) = quantile(abs(postMCMC),.975);  % upper bound on 2 sided CS
            cpostlow(j,i) = quantile(abs(postMCMC),.0275); % lower bound on 2 sided CS
            clear logpdf postMCMC
        end
        % Calculate quantiles for the bootstrap
        boot = g(bsxfun(@plus,Z'/(n^.5), reshape(thetaML(j,:), [mcov,1])));
        cboothigh(j, :) = reshape(quantile(boot, .975,2),[1,mcov]);
        cbootlow(j, :) = reshape(quantile(boot, .0275,2),[1,mcov]);
    end

    coveragepost = mean(bsxfun(@le,cpostlow,g(theta)).*...
        bsxfun(@le,g(theta),cposthigh),2);
    coverageboot = mean(bsxfun(@le,cbootlow,g(theta)).*...
        bsxfun(@le,g(theta),cboothigh),2);

    % 2.3.4.) Report the coverage of this procedure over theta at deltstar
    figure(ifig);
    plot(theta,coveragepost,'-.red'); hold on
    plot(theta,coverageboot,'-.blue'); hold on
    axis([-2 2 0 1.2])
    legend({'95% Credible Set based on the posterior quantiles','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11);
    legend('boxoff')
    xlabel('\theta')
    ylabel('Coverage')
    title(['Coverage Probability of 95\% CS and Bootstrap CI for $|\theta|$ with Prior $\theta \sim \gamma(2,2)-3$'],'Interpreter','latex');
    hold off
    saveTightFigure(figure(ifig), strcat('abs2_','n=',num2str(n),'Posterior95csGamma_Coverage.pdf'));
    ifig = ifig + 1;
end

%% ------------------------------------------------------------
% 3. Figures: Beta Prior
% ------------------------------------------------------------
disp('Gamma Prior: Theta 2*(width+0.5)(Z - 0.5), Z ~ Beta(alphap,betap)')

% 3.1. Set up functions:
% ------------------------------------------------------------
% %prior Theta = 2*(width+0.5)(Z - 0.5), Z ~ Beta(alphap,betap)
% % Note that the algorithm can be adapted to a different prior by just
% % changing definiton of the prior below and the function logpdfaux
% %parameters
alphap = 2;
betap = 2;
logpdfaux = @(x,data, sig0, alphap, betap) (-normlike([x,sig0],data) + log(betapdf((x/(2*(width + 0.5))+0.5),alphap, betap)));
logproppdf = @(x,y) -.5*((x-y).^2);
proprnd= @(x) x+randn;

% % 3.2a. Calculate quantiles for a set vector of data realisation:
% % ------------------------------------------------------------
if quantbeta_a == 1
    % 3.2.1.) Set up the containers
    thetadata=zeros(m,1);
    qgthetahigh_2S=zeros(m,1); % Container for upper bound on 2 sided CS
    qgthetalow_2S=zeros(m,1); % Container for lower bound on 2 sided CS

    for i=1:m
        i
        logpdf=@(x) logpdfaux(x,dataMC(i,:),sig0, alphap, betap);
        postMCMC = mhsample (0, I, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',100);
        % Generate I posterior draws from theta
        thetadata(i,:)=(n^.5)*mean(dataMC(i,:));            %thetahat
        % Calculate the quantiles
        qgthetahigh_2S(i,:)=quantile(abs(postMCMC),.975); % upper bound on 2 sided CS
        qgthetalow_2S(i,:)=quantile(abs(postMCMC),.025); % lower bound on 2 sided CS
        clear logpdf postMCMC
    end

    % 3.2.2.) Bootstrap and Asymptotic Approximation (Two Sided)
    figure(ifig)
    step=sthetaML(1,end)-sthetaML(1,end-1);
    h1=scatter(thetadata,qgthetahigh_2S,1.5,'red'); hold on
    scatter(thetadata,qgthetalow_2S,1.5,'red'); hold on
    for ix=1:size(sthetaML,2)
        h2=plot([sthetaML(1,ix)  sthetaML(1,ix)],[Approx_quant_2S(2,ix) Approx_quant_2S(1,ix)],':blue'); hold on
        plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(2,ix) Approx_quant_2S(2,ix)],'blue'); hold on
        plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(1,ix) Approx_quant_2S(1,ix)],'blue'); hold on
    end
    axis([-3.5 3.5 0 1])
    legend([h1,h2],{'End Points of 95% CS from posterior quantiles (MCMC)','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11)
    legend('boxoff')
    box on
    hold off
    xlabel('$Z_n = \sqrt{n} \hat{\theta}_n$','Interpreter','latex');
    title(['Credible Sets and Bootstrap CIs for $|\theta|$ with Prior $\theta \sim (\beta(2,2)-0.5)\times 5$'],'Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('abs1_','n=',num2str(n),'Posterior95csBeta_2Side_Rand.pdf'));
    ifig = ifig + 1;

end

% 3.2b. Calculate quantiles for a set vector of data realisation (using grid of MLE):
% ------------------------------------------------------------
normlike_explicit = @(x, thetaML, sig0) ( normpdf(x, thetaML, sig0/sqrt(n)) );
logpdfaux_explicit = @(x,thetaML, sig0, alphap, betap) (log(normlike_explicit(x, thetaML, sig0)) + log(betapdf((x/(2*(width + 0.5))+0.5),alphap, betap)));

if quantbeta_b == 1
    % 2.2b.1.) Set up the containers
    qgthetahigh_2S=zeros(lenML,1); % Container for upper bound on 2 sided CS
    qgthetalow_2S=zeros(lenML,1); % Container for lower bound on 2 sided CS

    for i=1:lenML
        i
        % Generate I posterior draws from posterior distribution
        logpdf=@(x) logpdfaux_explicit(x,sthetaML(i)/sqrt(n),sig0, alphap, betap);
        postMCMC = mhsample (0, IG, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',100);
        
        % Calculate the quantiles
        qgthetahigh_2S(i,:)=quantile(abs(postMCMC),.975); % upper bound on 2 sided CS
        qgthetalow_2S(i,:)=quantile(abs(postMCMC),.025); % lower bound on 2 sided CS
        clear logpdf postMCMC
    end

    % 2.2b.2.) Bootstrap and Asymptotic Approximation (Two Sided)
    figure(ifig)
    step=sthetaML(1,end)-sthetaML(1,end-1);
    for ix=1:size(sthetaML,2)
        h1=plot([sthetaML(1,ix)  sthetaML(1,ix)],[qgthetalow_2S(ix) qgthetahigh_2S(ix)],':red'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[qgthetahigh_2S(ix) qgthetahigh_2S(ix)],'red'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[qgthetalow_2S(ix) qgthetalow_2S(ix)],'red'); hold on
        
    	h2=plot([sthetaML(1,ix)  sthetaML(1,ix)],[Approx_quant_2S(2,ix) Approx_quant_2S(1,ix)],':blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(2,ix) Approx_quant_2S(2,ix)],'blue'); hold on
    	plot([sthetaML(1,ix)-(step/4)  sthetaML(1,ix)+(step/4)],[Approx_quant_2S(1,ix) Approx_quant_2S(1,ix)],'blue'); hold on
    end
    axis([-3.5 3.5 0 1])
    legend([h1,h2],{'End Points of 95% CS from posterior quantiles (MCMC)','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11)
    legend('boxoff')
    box on
    hold off
    xlabel('$Z_n = \sqrt{n} \hat{\theta}_n$','Interpreter','latex');
    title(['Credible Sets and Bootstrap CIs for $|\theta|$ with Prior $\theta \sim (\beta(2,2)-0.5)\times 5$'],'Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('abs1_','n=',num2str(n),'Posterior95csBeta_2Side.pdf'));
    ifig = ifig + 1;
end

% 3.3. Calculate covegage for a set vector of data realisation:
% ------------------------------------------------------------
if covbeta == 1
    % 3.3.1.) Calculate data realisations
    database = randn(mcov,n);
    thetaML = bsxfun(@plus,(1/(n^.5))*randn(1,mcov),theta); % MLE from data realisations (This is J x I)

    % 3.3.2.) Set up containers for quantiles
    cposthigh = zeros(J, mcov);
    cpostlow = zeros(J, mcov);
    cboothigh = zeros(J, mcov);
    cbootlow = zeros(J, mcov);

    % 3.3.3.) Calculate quantiles
    for j = 1:J
        j
        dataMCtheta = database + theta(j);
        for i=1:mcov
            % Generate I posterior draws from posterior distribution
            logpdf = @(x) logpdfaux(x,dataMCtheta(i,:),sig0, alphap, betap);
            postMCMC = mhsample(theta(j), Icov, 'logpdf', logpdf, 'logproppdf', logproppdf, 'proprnd', proprnd,'symmetric',1,'burnin',50);
            % Calcualte the quantiles
            cposthigh(j,i) = quantile(abs(postMCMC),.975);  % upper bound on 2 sided CS
            cpostlow(j,i) = quantile(abs(postMCMC),.0275); % lower bound on 2 sided CS
            clear logpdf postMCMC
        end
        % Calculate quantiles for the bootstrap
        boot = g(bsxfun(@plus,Z'/(n^.5), reshape(thetaML(j,:), [mcov,1])));
        cboothigh(j, :) = reshape(quantile(boot, .975,2),[1,mcov]);
        cbootlow(j, :) = reshape(quantile(boot, .0275,2),[1,mcov]);
    end

    coveragepost = mean(bsxfun(@le,cpostlow,g(theta)).*...
        bsxfun(@le,g(theta),cposthigh),2);
    coverageboot = mean(bsxfun(@le,cbootlow,g(theta)).*...
        bsxfun(@le,g(theta),cboothigh),2);

    % 3.3.4.) Report the coverage of this procedure over theta at deltstar
    figure(ifig);
    plot(theta,coveragepost,'-.red'); hold on
    plot(theta,coverageboot,'-.blue'); hold on
    axis([-2 2 0 1.2])
    legend({'95% Credible Set based on the posterior quantiles','95% Confidence Set based on the parametric Bootstrap'},'Location','Northwest','FontSize',11);
    legend('boxoff')
    xlabel('\theta')
    ylabel('Coverage')
    title(['Coverage Probability of 95\% CS and Bootstrap CI for $|\theta|$ with Prior $\theta \sim (\beta(2,2)-0.5)\times 5$'],'Interpreter','latex');
    hold off
    saveTightFigure(figure(ifig), strcat('abs2_','n=',num2str(n),'Posterior95csBeta_Coverage.pdf'));
    ifig = ifig + 1;
end
