function s=merge_em_spectrum(s1,s2)

% merge 2 spectra that have different x values
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2017

plot_spectra=0;

if isempty(s1); s=s2; return; end
if isempty(s2); s=s1; return; end

x1=s1(:,1);
y1=s1(:,2);
xs1=circshift(x1,-1)-x1;
d1=min(xs1(1:end-1));

x2=s2(:,1);
y2=s2(:,2);
xs2=circshift(x2,-1)-x2;
d2=min(xs2(1:end-1));

[x,sorted_x]=sort([x1;x2]);
tmp=[y1;y2];
y=tmp(sorted_x);

d=d1<=d2;
if d==1 % s1 has the best resolution
    for i=1:length(sorted_x)
        if sorted_x(i)>length(x1) % dealing with an element of s2
            if sorted_x(i)>1 && sorted_x(i)<length(x)
                [t,w_t]=min([(x(i)-x(i-1)),(x(i+1)-x(i))]);
                if t<d1
                    if w_t==1
                        x(i)=-1;
                        y(i-1)=y(i)+y(i-1);
                    else
                        x(i)=-1;
                        y(i+1)=y(i)+y(i+1);
                    end
                end
            end
        end
    end
else % s2 has the best resolution
    for i=1:length(sorted_x)
        if sorted_x(i)<=length(x1) % dealing with an element of s1
            if sorted_x(i)>1 && sorted_x(i)<length(x)
                [t,w_t]=min([(x(sorted_x(i))-x(sorted_x(i)-1)),(x(sorted_x(i)+1)-x(sorted_x(i)))]);
                if t<d1
                    if w_t==1
                        x(sorted_x(i))=-1;
                        y(sorted_x(i)-1)=y(sorted_x(i))+y(sorted_x(i)-1);
                    else
                        x(sorted_x(i))=-1;
                        y(sorted_x(i)+1)=y(sorted_x(i))+y(sorted_x(i)+1);
                    end
                end
            end
        end
    end
end

y(x==-1)=[];
x(x==-1)=[];
s=[x,y];

if plot_spectra==1
    subplot(2,1,1)
    plot(x,y)
    subplot(2,1,2)
    plot(x1,y1); hold on; plot(x2,y2); hold off
end

end

