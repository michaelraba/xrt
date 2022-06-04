% for radial correlation - classic pod
function plotSkmr(plotObject,array2,isGraph)
[ntimesteps rMin ,rMax, ss ,ncs ,plotOn ,azimuthalSet, azimuthalSetSize ,printStatus, lags,csSet,timeSet]=constants();
if isGraph=="graph"
  f=figure('Renderer', 'painters', 'Position', [10 10 1900 900])
  cMaxx=3
  A=linspace(0,1,540)
  %A=linspace(0,1,1079)
  xlabel('radius $\frac{r}{R}$','interpreter','latex')
  ylabel("$S_{ii}(k,m;r,r')$",'interpreter','latex')  
  %xlabel('radius r');
  %ylabel('S');
  %figure;
  hold on;
  cou = 1
  for c=1:cMaxx
      subplot(3,1,c);

      %for t=1:ntimesteps
      for m=2:18
      
      labelStr = ['(m,k)=(', num2str(azimuthalSet(m)),',',num2str(c),')'];
      %pp=plot(A,real(plotObject(c).circle(m).dat )/ntimesteps,"DisplayName",labelStr);
      pp=plot(real(plotObject(c).circle(m).dat(1:540)/ntimesteps ),"DisplayName",labelStr);


          tiSt=['$\Big\langle\Phi^{(' num2str(c ) ')}_{x}(m;r)\Big\rangle_k$'];
          title(tiSt, 'FontName','capitana','FontSize',12,'interpreter','latex')
      if c==1      
      legend();
      end
      %avgTimeEnd(1).circle(1).dat  
      %pp=plot(A,real(plotObject(c).circle(m).dat /ntimesteps),"DisplayName",labelStr)

      
      
      %pp=plot(A,real(plotObject(c).t(t).azimuthal(m).dat((end-1)/2:end) )/ntimesteps,"DisplayName",labelStr)
                    %SrrPrForFourierPost(c).t(t).azimuthal(m).dat
      %if mod(c,4)==0
      %  pp.LineStyle='-';
      %elseif mod(c,4)==1
      %  pp.LineStyle='--';
      %elseif mod(c,4)==2
      %  pp.LineStyle=':';
      %elseif mod(c,4)==3
      %  pp.LineStyle='-.';
      %end
      hold on;
      %pause(1);
      end % m
      cou = cou + 1;
  end %c
     % end %t
      %legend( );
      %saveas(gcf,'Sij.png','res')
      %sprintf('%s','pause');
elseif isGraph=="graphPause"
  hold on;
    for m=2:38
    for c=1:ncs
    plot(real(plotObject(m).crosssection(c).dat(end/2:end) )/ntimesteps)
    hold on;
    pause(1);
    end
    end
    %sprintf('%s','pause');
      pause(1)
%else
%yy=csSet(2);
kStr= [', $k \in [1,$' num2str(csSet(2)) '$,\ldots,$' num2str(csSet(ncs))  '] (k-modes are averaged).'];
tStr= [ ', $t\in[1,' num2str(timeSet(2)) ',' num2str(timeSet(3)) ',' num2str(timeSet(4)) ',\ldots,' num2str(timeSet(ntimesteps))  ',...,999]$' ]
aziStr= [' , $m\in [1,$' num2str(azimuthalSet(2)) '$,\ldots,$' num2str(azimuthalSet(azimuthalSetSize)) ']'];
titleStrr=['Classic POD modes $\Big\langle\Phi_{x}(m;r)\Big\rangle_k$ '   tStr kStr]
sgtitle(titleStrr,'FontName','capitana','FontSize',12,'interpreter','latex')
 


%% graph timeAvg

elseif isGraph=="timeAvg"
  f=figure('Renderer', 'painters', 'Position', [10 10 1900 900])
  cMaxx=3
  %A=linspace(0,1,540)
  A=linspace(0,1,1079)
  maxT = 1;
  xlabel('radius $\frac{r}{R}$','interpreter','latex')
  ylabel("$S_{ii}(k,m;r,r')$",'interpreter','latex')  
  %xlabel('radius r');
  %ylabel('S');
  %figure;
  hold on;
  cou = 1
  for c=1:1
      subplot(1,2,c);
      %for t=1:ntimesteps
      for m=1:18
      labelStr = ['(m,k)=(', num2str(azimuthalSet(m)),',',num2str(c),')'];
      %pp=plot(A,real(plotObject(c).circle(m).dat((end-1)/2:end) )/ntimesteps,"DisplayName",labelStr);
      %pp=plot(A,real(plotObject(c).circle(m).dat )/ntimesteps,"DisplayName",labelStr);
      %pp=plot(A,real(plotObject.circle(m).dat )/ntimesteps,"DisplayName",labelStr);
      pp=plot(real(plotObject(c).circle(m).dat(1:540) )/ntimesteps,"DisplayName",labelStr);

          tiSt=["$R(m;k;r,r')$ Correlation With fft-x and time-averaging"];
          title(tiSt, 'FontName','capitana','FontSize',12,'interpreter','latex')
      if c==1      
      legend();
      end
      hold on;
      end % m
      cou = cou + 1;


subplot(1,2,2);
      for t=1:maxT
      
      for m=1:18
      labelStr = ['(m,k)=(', num2str(azimuthalSet(m)),',',num2str(c),')'];
      %pp=plot(A,real(plotObject(c).circle(m).dat )/ntimesteps,"DisplayName",labelStr);
      plot(real(array2(t).circle(m).dat(540:end) ),"DisplayName",labelStr);
      %xcorrDone(9).circle(18).dat  
          tiSt=["$R(m;x;r,r')$ Correlation without fft-x; shows a particular snapshot t=1."];
          title(tiSt, 'FontName','capitana','FontSize',12,'interpreter','latex')
      if c==1      
      legend();
      end
      hold on;
      end % m
      cou = cou + 1;
      end 


  end %c
aziStr= ["Correlation $R(k;m;r,r')$ (time-averaged)" ];
 % titleStrr=['Classic POD modes $\Big\langle\Phi_{x}(m;r)\Big\rangle_k$ for (tTot,xTot)=('  aziStr  kStr]
  sgtitle(aziStr,'FontName','capitana','FontSize',12,'interpreter','latex')
end % if

end % f
