function psiplot(psi)
global x y N a b 
pcolor(x,y,reshape(abs(psi).^2,N,N));
%pcolor(x,y,reshape(angle(psi),N,N));
colormap('parula');
axis equal
xlim([a,b]);ylim([a,b]);
shading interp;
colorbar;
drawnow;
end

