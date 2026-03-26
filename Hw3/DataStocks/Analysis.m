% Preliminary code to transform data into matrix after loading
% 
D(:,:,1)=AAPL(2:end,:);
D(:,:,2)=ADBE(2:end,:);
D(:,:,3)=AMZN(2:end,:);
D(:,:,4)=AVGO(2:end,:);
D(:,:,5)=FB(2:end,:);
D(:,:,6)=GOOGL(2:end,:);
D(:,:,7)=MSFT(2:end,:);
D(:,:,8)=QCOM(2:end,:);
D(:,:,9)=TWTR(2:end,:);
D(:,:,10)=AMD(2:end,:);
date=[1:length(D(:,1,1))]';
T=['Open  ';'High  '; 'Low   '; 'Close '; 'Volume'];
Stock=['AAPL';'ADBE';'AMZN';'AVGO';'FB  ';'GOOG';'MSFT';'QCOM';'TWTR';'AMD '];

figure(1)
d=(D(:,2,:)-D(:,3,:))./D(:,1,:);
d=squeeze(d);
for k=1:10,
    subplot(10,1,k)
    plot(date,d(:,k))
    title(['Intra-day dispersion ',Stock(k,:),'  mean=',num2str(100*mean(d(:,k)))])
end

figure(2)
d=(D(:,2,:)-D(:,3,:))./D(:,1,:);
d=squeeze(d);

for k=1:10,
    subplot(5,2,k)
    plot(squeeze(D(:,5,k)),100*d(:,k),'.')
    title(['Dispersion vs volume ',Stock(k,:),'  mean=',num2str(100*mean(d(:,k)))])
    xlabel('Volume')
    ylabel('Disp. %')
end

