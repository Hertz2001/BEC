function [ngrad, energy, mu]=becfff(psi0)
global bet omg N area h v2d x y

%input is a vector of size (2*N*N,1)
% psireal=reshape(psi0(1:(N-2)*(N-2)),N-2,N-2);
% psiimag=reshape(psi0((N-2)*(N-2)+1:2*(N-2)*(N-2)),N-2,N-2);
% psire=zeros(N);psiim=zeros(N);
% psire(2:N-1,2:N-1)=psireal;psiim(2:N-1,2:N-1)=psiimag;

psire=reshape(psi0(1:N*N),N,N);
psiim=reshape(psi0(N*N+1:2*N*N),N,N);


%lapla=([psi(:,N) psi(:,1:N-1)]+[psi(:,2:N) psi(:,1)]+[psi(N,:);psi(1:N-1,:)]+[psi(2:N,:);psi(1,:)]-4*psi)/(h*h);
laplare=([psire(:,N) psire(:,1:N-1)]+[psire(:,2:N) psire(:,1)]+[psire(N,:);psire(1:N-1,:)]+[psire(2:N,:);psire(1,:)]-4*psire)/(h*h);
laplaim=([psiim(:,N) psiim(:,1:N-1)]+[psiim(:,2:N) psiim(:,1)]+[psiim(N,:);psiim(1:N-1,:)]+[psiim(2:N,:);psiim(1,:)]-4*psiim)/(h*h);

%gradx=[psi(2,:)-psi(N,:);psi(3:N,:)-psi(1:N-2,:);psi(1,:)-psi(N-1,:)]/(2*h);
gradxre=[psire(2,:)-psire(N,:);psire(3:N,:)-psire(1:N-2,:);psire(1,:)-psire(N-1,:)]/(2*h);
gradxim=[psiim(2,:)-psiim(N,:);psiim(3:N,:)-psiim(1:N-2,:);psiim(1,:)-psiim(N-1,:)]/(2*h);

ydxre=y .*gradxim;
ydxim=y .*gradxre;
%grady=[psi(:,2)-psi(:,N) psi(:,3:N)-psi(:,1:N-2) psi(:,1)-psi(:,N-1)]/(2*h);
gradyre=[psire(:,2)-psire(:,N) psire(:,3:N)-psire(:,1:N-2) psire(:,1)-psire(:,N-1)]/(2*h);
gradyim=[psiim(:,2)-psiim(:,N) psiim(:,3:N)-psiim(:,1:N-2) psiim(:,1)-psiim(:,N-1)]/(2*h);

xdyre=x' .*gradyim;
xdyim=x' .*gradyre; 
rotatre=(xdyre-ydxre);
rotatim=ydxim-xdyim;

potenre=v2d.*psire;
potenim=v2d.*psiim;

psilengthsquare=psire.^2+psiim.^2;
cubicre=psilengthsquare.*psire;
cubicim=psilengthsquare.*psiim;

%ngrad=lapla/2+omg*rotat-poten-bet*cubic;
ngradre=laplare/2+omg*rotatre-potenre-bet*cubicre;
%ngradre=-(h*h)*laplare-2*h*h*omg*rotatre+2*h*h*potenre+2*h*h*bet*cubicre
ngradim=laplaim/2+omg*rotatim-potenim-bet*cubicim;
%ngradim=-(h*h)*laplaim-2*h*h*omg*rotatim+2*h*h*potenim+2*h*h*bet*cubicim

if nargout>1
    int_g=norm(gradxre,'fro')^2+norm(gradyre,'fro')^2+norm(gradxim,'fro')^2+norm(gradyim,'fro')^2;
    int_r=(sum(sum(-rotatre.*psire-rotatim.*psiim)));
    int_p=(sum(sum(v2d.*psilengthsquare)));
    int_n=(sum(sum(psilengthsquare.^2)));
    energy=(int_g/2+omg*int_r+int_p+bet/2*int_n)*area;
   
    mu=energy+bet/2*int_n*area;

%     ngrad=ngrad+mu*psi;%/normpsi
end
%ngradre=ngradre(2:N-1,2:N-1);ngradim=ngradim(2:N-1,2:N-1);
ngradre=ngradre(:);ngradim=ngradim(:);
ngrad=[ngradre;ngradim];
end
