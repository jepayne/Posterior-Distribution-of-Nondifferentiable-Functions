% ----------------------------------------------------------------------------
% Max Value Example
% ----------------------------------------------------------------------------
% Description:  This file comptutes the plots for Figures 1 and 2 in the
%   online appendix to Kitagawa, Montiel-Olea, Payne. "Posterior Distribution of
%   Non-Differential Functions". (2016).
% ----------------------------------------------------------------------------
% Authors:      Toru Kitagawa (t.kitagawa@ucl.ac.uk),
%               Jose-Luis Montiel-Olea (montiel.olea@nyu.edu), and
%               Jonathan Payne (jep459@nyu.edu)
% ----------------------------------------------------------------------------


%% ----------------------------------------------------------------------------
% 0. Set up
% ----------------------------------------------------------------------------

clear; clc;
rng('default') % Set random seed to default

% Define transformations
g=@(theta1,theta2)(max(theta1,theta2));
gprime=@(h1,h2,theta1,theta2)((theta1>theta2)*h1+(theta1<theta2)*h2+(theta1==theta2)*max(h1,h2));

% Posterior precision
lambda2=1;

% Parameters from statistical model
sigma = [1,0;0,1];

% Loop parameters
n = 100; % Sample size
nZ = 100; % Number of draws to estimate parametric bootstrap
I = 1000; % Number of samples [Set to 10000 in paper]
IZn = 1000; % Number of draws from posterior [Set to 10000 in paper]
J = 101; % Number of points in grid for theta
width = 2;

ifig = 1;

% Generate random variables for sampling from the posterior or bootstrap
Z = mvnrnd(zeros(2,1),eye(2), nZ)';

% Adjustment parameter for Hong method
delta=1/4;

% Indicator variables (to indicate which sections of the code are run)
covtheta12equal = 1;
covHong = 1;
credHong = 1;

%% ------------------------------------------------------------
% Figure 1: Coverage probabilities for:
%   \theta_1 = \theta_2
% ------------------------------------------------------------
if covtheta12equal == 1

    for nt = [100, 200, 300, 1000]
        % 1.1. Set up containers for quantiles
        cposthigh = zeros(J, I); % Container for upper quantiles of posterior
        cpostlow = zeros(J, I); % Container for lower quantiles of posterior
        cboothigh = zeros(J, I); % Container for upper quantiles of bootstrap
        cbootlow = zeros(J, I); % Container for lower quantiles of bootstrap

        % 1.2. Generate data sequence
        theta = [linspace(-width,width,J)', linspace(-width,width,J)']; % Grid for the parameter space
        sigma=[1,0;0,1];
        thetaML = zeros(2,I,J);
        Zstar = mvnrnd(zeros(2,1), sigma, I)';

        % 1.3. Calculate quantiles
        for j = 1:J
            j
            % Generate data
            thetaML(:,:,j) = bsxfun(@plus, nt^(-0.5)*Zstar, theta(j,:)');
            % Calculate quantiles for the posterior
            post = g(bsxfun(@plus,Z(1,:)/((nt+lambda2)^.5), (nt/(nt+lambda2))*reshape(thetaML(1,:,j), [I,1])),...
                bsxfun(@plus,Z(2,:)/((nt+lambda2)^.5), (nt/(nt+lambda2))*reshape(thetaML(2,:,j), [I,1])));
            cposthigh(j, :) = reshape(quantile(post, .975,2),[1,I]);
            cpostlow(j, :) = reshape(quantile(post, .0275,2),[1,I]);
            % Calculate quantiles for the bootstrap
            boot = g(bsxfun(@plus,Z(1,:)/((nt)^.5), reshape(thetaML(1,:,j), [I,1])),...
                bsxfun(@plus,Z(2,:)/((nt)^.5), reshape(thetaML(2,:,j), [I,1])));
            cboothigh(j, :) = reshape(quantile(boot, .975,2),[1,I]);
            cbootlow(j, :) = reshape(quantile(boot, .0275,2),[1,I]);
        end
        coveragepost = mean(bsxfun(@le,cpostlow,g(theta(:,1),theta(:,2))).*...
            bsxfun(@le,g(theta(:,1),theta(:,2)),cposthigh),2);
        coverageboot = mean(bsxfun(@le,cbootlow,g(theta(:,1),theta(:,2))).*...
            bsxfun(@le,g(theta(:,1),theta(:,2)),cboothigh),2);

        % 1.4. Report the coverage
        figure(ifig);
        plot(theta(:,1),coveragepost,'-.red'); hold on
        plot(theta(:,1),coverageboot,'-.blue'); hold on
        axis([-2 2 0 1])
        legend({'95% Credible Set based on the posterior quantiles','95% Confidence Set based on the parametric Bootstrap'},'Location','Southwest','FontSize',11);
        legend('boxoff')
        xlabel('\theta')
        ylabel('Coverage')
        hold off
        saveTightFigure(figure(ifig), strcat('max4_','lambda2=',num2str(lambda2),'_n=',num2str(nt),'_M=',num2str(theta(end,1)),'_CoverageThetaEqual.pdf'));
        ifig = ifig + 1;
    end
end

%% ------------------------------------------------------------
% Figures 2: Hong Coverage
% Heat map: \theta_1 \in (-width, width) and \theta_2 \in (-width, width)
% ------------------------------------------------------------

if covHong == 1
    % 2.1. Set up containers for quantiles
    chonghigh = zeros(J, I); % Container for upper quantiles of posterior
    chonglow = zeros(J, I); % Container for lower quantiles of posterior
    CovMatrix = zeros(J, J);

    % 2.2. Generate data sequence
    theta1 = linspace(-width,width,J)';
    theta2 = linspace(-width,width,J)';
    sigma=[1,0;0,1];
    thetaML = zeros(2,I,J);
    thetaMLij = zeros(2,I);
    Zstar = mvnrnd(zeros(2,1), sigma, I)';

    sigma=[1,0;0,1];
    Znstar=mvnrnd(zeros(2,1), sigma, IZn);

    delta=1/4;

    % 2.3. Calculate quantiles
    for j1 = 1:J
        j1
        for j2 = 1:J
            % Generate data given the theta = [theta1(j1); theta2(j2)]
            thetaMLij = bsxfun(@plus, n^(-0.5)*Zstar, [theta1(j1); theta2(j2)]);
            % Calculate Hong approximate random variable
            Honginput1 = bsxfun(@plus, thetaMLij(1,:),Znstar(:,1)/(n^(0.5-delta))); % size = (IZn x I)
            Honginput2 = bsxfun(@plus, thetaMLij(2,:),Znstar(:,2)/(n^(0.5-delta))); % size = (IZn x I)
            Hongapprox = n^(0.5-delta)*bsxfun(@plus, g(Honginput1,Honginput2),-g(thetaMLij(1,:),thetaMLij(2,:)));
            % Calculate quantiles for the hong approach
            gthetahat = g(thetaMLij(1,:), thetaMLij(2,:));
            cvs=bsxfun(@plus,-n^(-0.5)*quantile(Hongapprox,[.025,.975],1), gthetahat);
            CovMatrix(j1,j2) = mean(bsxfun(@le,cvs(2,:),g(theta1(j1),theta1(j2))).*...
                bsxfun(@le,g(theta1(j1),theta1(j2)),cvs(1,:)),2);
        end
    end

    % 2.4. Report the coverage
    figure(ifig);
    colormap('hot');
    scale = [0.6 1];
    imagesc(theta1,theta2,CovMatrix, scale);
    set(gca,'YDir','normal')
    colorbar;
    xlabel('$\theta_1$','Interpreter','latex');
    ylabel('$\theta_2$','Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('max5_','lambda2=',num2str(lambda2),'_n=',num2str(n),'_CovHongHeat.pdf'));
    ifig = ifig + 1;

end

%% ------------------------------------------------------------
% Figures 3: Hong Credibility (Using Posterior Distribution)
% Heat map: \theta_1 \in (-width, width) and \theta_2 \in (-width, width)
% ------------------------------------------------------------
if credHong == 1
    thetahat1 = linspace(-2,2,J)'; % Possible values for thetahat1
    thetahat2 = linspace(-2,2,J)'; % Possible values for thetahat2
    CredMatrix = zeros(J,J); % Will hold credibility values

    for j2 = 1:J
        j2
        % Define thetahat
        thetahat = [thetahat1, thetahat2(j2)*ones(1,J)']';

        % Calculate posterior quantiles
        Znstar = mvnrnd(zeros(2,1), sigma, I); % Draws of Znstar for calculating the Hong quantiles
        Honginput1 = bsxfun(@plus, (n^(delta)*lambda2/(n+lambda2) + 1)*thetahat(1,:),Znstar(:,1)/(n^(-delta)*(n+lambda2)^(0.5)));
        Honginput2 = bsxfun(@plus, (n^(delta)*lambda2/(n+lambda2) + 1)*thetahat(2,:),Znstar(:,2)/(n^(-delta)*(n+lambda2)^(0.5)));
        Hongapprox = n^(0.5-delta)*bsxfun(@plus, g(Honginput1,Honginput2),-g(thetahat(1,:),thetahat(2,:)));
        gthetahat = g(thetahat(1,:), thetahat(2,:));
        cvs = bsxfun(@plus,-n^(-0.5)*quantile(Hongapprox,[.025,.975],1), gthetahat);

        % Generate draws from the posterior
        ZP = mvnrnd(zeros(2,1), sigma, IZn); % Draws to build the posterior
%         ZP = Znstar;
        input1 = bsxfun(@plus, thetahat(1,:)*n/(n+lambda2), ZP(:,1)*(n+lambda2)^(-0.5));
        input2 = bsxfun(@plus, thetahat(2,:)*n/(n+lambda2), ZP(:,2)*(n+lambda2)^(-0.5));
        gpost = g(input1,input2);

        CredMatrix(:,j2)=mean(bsxfun(@le,gpost,cvs(1,:)).*bsxfun(@ge,gpost,cvs(2,:)))';
    end

    figure(ifig)
    colormap('hot');
    scale = [0.6 1];
    imagesc(thetahat1, thetahat2, CredMatrix, scale);
    set(gca,'YDir','normal')
    colorbar;
    xlabel('$\hat{\theta}_1$','Interpreter','latex');
    ylabel('$\hat{\theta}_2$','Interpreter','latex');
    saveTightFigure(figure(ifig), strcat('max6_','lambda2=',num2str(lambda2),'_n=',num2str(n),'_M=',num2str(thetahat(end,1)),'_CredThetaHeatboth.pdf'));
    ifig = ifig + 1;

end
