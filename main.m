close all;
clear all;

headers();

% =========================================================================

omega = 3;
tau = 1;
stepTab = [0.5 0.25 0.125 0.0625];
%stepTab = [0.5 0.25 0.125 0.0625 0.0312];

errorL2CG1   = zeros(1,size(stepTab,2));
errorH1CG1   = zeros(1,size(stepTab,2));
condACG1     = zeros(1,size(stepTab,2));
errorL2CG2   = zeros(1,size(stepTab,2));
errorH1CG2   = zeros(1,size(stepTab,2));
condACG2     = zeros(1,size(stepTab,2));
errorL2HDG1 = zeros(1,size(stepTab,2));
errorH1HDG1 = zeros(1,size(stepTab,2));
condHDG1    = zeros(1,size(stepTab,2));
condHDG1red = zeros(1,size(stepTab,2));
errorL2HDG2 = zeros(1,size(stepTab,2));
errorH1HDG2 = zeros(1,size(stepTab,2));
condHDG2    = zeros(1,size(stepTab,2));
condHDG2red = zeros(1,size(stepTab,2));
errorL2DG = zeros(1,size(stepTab,2));
errorH1DG = zeros(1,size(stepTab,2));
condDG    = zeros(1,size(stepTab,2));
condDGred = zeros(1,size(stepTab,2));

for k=1:size(stepTab,2)
    
    step = stepTab(k);
    system(['gmsh -2 mesh.geo -clmax ' num2str(step) ' -clmin ' num2str(step)]);
    readMesh('mesh.msh');
    
    buildMatrixGloCG();
    [errorL2CG1(k), errorH1CG1(k), condACG1(k)] = computeSolNumCG1(@mySol2,omega);
    [errorL2CG2(k), errorH1CG2(k), condACG2(k)] = computeSolNumCG2(@mySol2,omega);
    
    buildMatrixGloHDG();
    [errorL2HDG1(k), errorH1HDG1(k), condHDG1(k), condHDG1red(k)] = computeSolNumHDG1(@mySol2,omega,tau,0);
    [errorL2HDG2(k), errorH1HDG2(k), condHDG2(k), condHDG2red(k)] = computeSolNumHDG1(@mySol2,omega,tau,1);
    
    %     buildMatrixGloDG();
    %     [errorL2DG(k), errorH1DG(k), condDG(k), condDGred(k)] = computeSolNumDG1(@mySol2,omega);
    
end

% global AA
% figure(3);
% spy(AA);

figure(1);
subplot(1,2,1);
loglog(1./stepTab,real(errorL2CG1),'k','LineWidth',2,'Marker','o','DisplayName','CG1 - L2 error');
hold on
loglog(1./stepTab,real(errorH1CG1),'k--','LineWidth',2,'Marker','o','DisplayName','CG1 - H1 error');
loglog(1./stepTab,real(errorL2CG2),'g','LineWidth',2,'Marker','o','DisplayName','CG2 - L2 error');
loglog(1./stepTab,real(errorH1CG2),'g--','LineWidth',2,'Marker','o','DisplayName','CG2 - H1 error');
loglog(1./stepTab,real(errorL2HDG1),'b','LineWidth',2,'Marker','*','DisplayName','HDG1 - L2 error');
loglog(1./stepTab,real(errorH1HDG1),'b--','LineWidth',2,'Marker','*','DisplayName','HDG1 - H1 error');
loglog(1./stepTab,real(errorL2HDG2),'r','LineWidth',2,'Marker','*','DisplayName','HDG2 - L2 error');
loglog(1./stepTab,real(errorH1HDG2),'r--','LineWidth',2,'Marker','*','DisplayName','HDG2 - H1 error');
loglog(1./stepTab,real(errorL2DG),'m','LineWidth',2,'Marker','*','DisplayName','DG - L2 error');
loglog(1./stepTab,real(errorH1DG),'m--','LineWidth',2,'Marker','*','DisplayName','DG - H1 error');
title('Omega 3 - Tau 1');
legend();
xlabel('1/h');
ylabel('Error');
subplot(1,2,2);
loglog(1./stepTab,condACG1,'k','LineWidth',2,'Marker','o','DisplayName','CG1');
hold on
loglog(1./stepTab,condACG2,'m','LineWidth',2,'Marker','o','DisplayName','CG2');
loglog(1./stepTab,condHDG1,'b','LineWidth',2,'Marker','*','DisplayName','HDG1');
loglog(1./stepTab,condHDG1red,'b--','LineWidth',2,'Marker','*','DisplayName','HDG1 (reduced)');
loglog(1./stepTab,condHDG2,'r','LineWidth',2,'Marker','*','DisplayName','HDG2');
loglog(1./stepTab,condHDG2red,'r--','LineWidth',2,'Marker','*','DisplayName','HDG2 (reduced)');
loglog(1./stepTab,condDG,'m','LineWidth',2,'Marker','*','DisplayName','DG');
loglog(1./stepTab,condDGred,'m--','LineWidth',2,'Marker','*','DisplayName','DG (reduced)');
title('Omega 3 - Tau 1');
legend();
xlabel('1/h');
ylabel('Conditionning');

% =========================================================================

% path1 = getenv('PATH');
% path1 = [path1 ':/Applications/Gmsh.app/Contents/MacOS/'];
% setenv('PATH', path1);
% step = 0.1;
% system(['gmsh -2 mesh.geo -clmax ' num2str(step) ' -clmin ' num2str(step)]);
%
% readMesh('mesh.msh');
% % vizuMesh();
%
% pbmFemCgBuild();
% omegaTab = 3; %0:0.05:10;
% errorL2cg = zeros(1,size(omegaTab,2));
% errorH1cg = zeros(1,size(omegaTab,2));
% condAcg = zeros(1,size(omegaTab,2));
% for k = 1:size(omegaTab,2)
%     omega = omegaTab(k);
%     [errorL2cg(k), errorH1cg(k), condAcg(k)] = pbmFemCgSolveSystem(@mySol2,omega);
% end
%
% % pbmFemHDgBuild();
% % omegaTab = 3; %0:0.05:10;
% % tau = 0.5;
% % errorL2 = zeros(1,size(omegaTab,2));
% % errorH1 = zeros(1,size(omegaTab,2));
% % condA = zeros(1,size(omegaTab,2));
% % for k = 1:size(omegaTab,2)
% %     omega = omegaTab(k);
% %     [errorL2(k), errorH1(k), condA(k)] = pbmFemHDgSolveSystem(@mySol2,omega,tau);
% % end
%
% pbmFemHDgBuild();
% omega = 3;
% tauTab = 0:0.05:2;
% errorL2hdg = zeros(1,size(tauTab,2));
% errorH1hdg = zeros(1,size(tauTab,2));
% condAhdg = zeros(1,size(tauTab,2));
% condAhdgEff = zeros(1,size(tauTab,2));
% for k = 1:size(tauTab,2)
%     tau = tauTab(k);
%     [errorL2hdg(k), errorH1hdg(k), condAhdg(k), condAhdgEff(k)] = pbmFemHDgSolveSystem(@mySol2,omega,tau);
% end
%
% % global AA
% % figure(3);
% % spy(AA);
%
% figure(1);
% subplot(1,2,1);
% semilogy(tauTab,real(errorL2hdg));
% hold on
% semilogy(tauTab,real(errorH1hdg));
% semilogy(tauTab,errorL2cg*ones(1,size(tauTab,2)))
% semilogy(tauTab,errorH1cg*ones(1,size(tauTab,2)))
% subplot(1,2,2);
% plot(tauTab,condAhdg);
% hold on
% plot(tauTab,condAhdgEff);
% plot(tauTab,condAcg*ones(1,size(tauTab,2)))

% =========================================================================

% figure(1);
% subplot(1,3,1);
% hold on;
% semilogy([1*pi/2 1*pi/2],[1e-15 1e5],'k:');
% semilogy([2*pi/2 2*pi/2],[1e-15 1e5],'k:');
% semilogy([3*pi/2 3*pi/2],[1e-15 1e5],'k:');
% semilogy([4*pi/2 4*pi/2],[1e-15 1e5],'k:');
% semilogy([5*pi/2 5*pi/2],[1e-15 1e5],'k:');
% semilogy([6*pi/2 6*pi/2],[1e-15 1e5],'k:');
% axis([0 10 1e-10 1e5]);
% xlabel('Wavenumber');
% ylabel('L2-error');
% subplot(1,3,2);
% hold on;
% semilogy([1*pi/2 1*pi/2],[1e-15 1e5],'k:');
% semilogy([2*pi/2 2*pi/2],[1e-15 1e5],'k:');
% semilogy([3*pi/2 3*pi/2],[1e-15 1e5],'k:');
% semilogy([4*pi/2 4*pi/2],[1e-15 1e5],'k:');
% semilogy([5*pi/2 5*pi/2],[1e-15 1e5],'k:');
% semilogy([6*pi/2 6*pi/2],[1e-15 1e5],'k:');
% axis([0 10 1e-10 1e5]);
% xlabel('Wavenumber');
% ylabel('H1-error');
% subplot(1,3,3);
% hold on;
% semilogy([1*pi/2 1*pi/2],[1 1e10],'k:');
% semilogy([2*pi/2 2*pi/2],[1 1e10],'k:');
% semilogy([3*pi/2 3*pi/2],[1 1e10],'k:');
% semilogy([4*pi/2 4*pi/2],[1 1e10],'k:');
% semilogy([5*pi/2 5*pi/2],[1 1e10],'k:');
% semilogy([6*pi/2 6*pi/2],[1 1e10],'k:');
% axis([0 10 1 1e10]);
% xlabel('Wavenumber');
% ylabel('Conditionning');
