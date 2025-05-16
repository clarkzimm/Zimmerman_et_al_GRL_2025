% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Computes analytic equilibrium hysteresis (DH) in hosing (H) for varying 
% gyre strengths (KN,KS), and boundary between monostable and bistable 
% (KN,KS) regimes. Called by Fig3_Kspace_DH_stats.m, 
% Fig4_Kspace_5_2x_diff.m and SI_FigS6_dethyst.m
%
% Dependencies:
%   - solve_H_3box.m OR solve_H_5box.m  (analytic equilibrium solutions of hosing 
%     H as a function of transport q) for 3-box or 5-box model
%   - parameters_Xbox_YY_ZCO2.m (parameter values for the X = [3 5]-box 
%       model as calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) 
%       with Z = [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%          
% Output:
%   - DH(KN,KS) (bistable solution width field as functions of KN, KS vectors)
%   - DHbound(KN) (parameterized boundary between mono and bistable 
%     regimes (i.e. DH(KN,KS)=0), given as KS0(KN))
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------
% set transport vector q for which equilibrium H is computed
qpos = linspace(0,20,1e3);  % q>0
qneg = linspace(-20,0,1e3); % q<0

% % uncomment if running here, leave commented if calling from
% % Fig3_Kspace_DH_stats.m or Fig4_Kspace_5_2x.m scripts
% % set model formulation (3 -> 3-box; 5 -> 5-box)
% box = 3;
% 
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM2-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;

%load parameters
params = sprintf('parameters_%dbox_%s_%dCO2',box,model,CO2);
run(params)

% load model solutions and parameters
solve_H = sprintf('solve_H_%dbox',box); 
eval(solve_H)

% compute DH for varying KN, KS
res = 1000; %length of KN,KS vectors; set to 1000 for production level figs; set to 100 for Figure S6b
KNvec = linspace(0,50,res); 
KSvec = linspace(0,50,res); 
DH = nan(length(KNvec),length(KSvec));
tic
for k = 1:length(KNvec)
    for l = 1:length(KSvec)
        
        KN = KNvec(k);
        KS = KSvec(l);

        Hp = Hpos(qpos,KN,KS); %Hpos is q>0 soln loaded from Hsol
        Hmax = Hp(Hp==max(Hp)); %Hmax corresponds to upper saddle node
        Hn = Hneg(qneg,KN,KS); %Hneg is q<0 soln loaded from Hsol
        Hmin = Hn(Hn==min(Hn)); %Hmin corresponds to lower saddle node

        %hysteresis as measured in distance between upper and lower saddle
        %nodes:
        DH(k,l) = Hmax-Hmin;
    end
end
toc
DH(DH<1e-10)=nan; %remove small numerical hysteresis

% compute boundary between monostable and bistable regime in KN,KS space
DHbound = nan(length(KNvec),1);
for i = 2:length(KNvec)
    DHrow = DH(i,:);
    if ~all(DHrow>0)
        DHbound(i) = KSvec(find(DHrow>0,1,'last'));
    end
end

%% save output
fname = sprintf('analytic_DH_%dbox_%s_%dCO2.mat',box,model,CO2);
save(fname,'KNvec','KSvec','DH','DHbound')
%% plot results (cf also Fig3_Kspace_DH_stats.m)
% figure(3); clf
% levs = 0:0.025:.4;
% cmap = colormap((brewermap(length(levs)-1,'YlGnBu')));
% contourf(KNvec,KSvec,DH',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
% hold on
% cc = colorbar;
% plot(KNvec,DHbound,'k--','linewidth',2,'HandleVisibility','off')
% clim([0 0.4])