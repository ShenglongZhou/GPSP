function RecoverShow(xo,x,pos)
    figure('Renderer', 'painters', 'Position', pos)
    axes('Position', [0.05 0.1 0.9 0.8] );
    stem(find(xo),xo(xo~=0),'bo-','MarkerSize',6, 'LineWidth',1),hold on
    stem(find(x),x(x~=0),'r*:', 'MarkerSize',4, 'LineWidth',1),hold on
    grid on, ymin= -0.1; ymax=0.2;
    xx  = [xo; x];
    if nnz(xx<0)>0, ymin= min(xx(xx<0))-0.1; end
    if nnz(xx>0)>0, ymax= max(xx(xx>0))+0.1; end   
    axis([1 length(x)  ymin  ymax])
    st1   = strcat('SNR=',num2str(-10*log10(norm(x-xo)^2),4));
    wrong = nnz(xo)-nnz(find(x~=0 & xo~=0)); 
    st2   = strcat(', s_*=',num2str(nnz(xo)),...
            ', Number of mis-supports =',num2str(wrong));         
    title(strcat(st1,st2))
    set(0,'DefaultAxesTitleFontWeight','normal');
    legend('Ground-Truth', 'Recovered')
end
