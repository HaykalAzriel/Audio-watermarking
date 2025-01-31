% function audiowatermarking1
% SWT-DCT-SVD-QIM
clear all
close all

%%%%%%%%%%%%%%%%%%
%Tidak terpakai (jangan diutak atik)
%%%%%%%%%%%%%%%%%%
M=64;N=64; %Resolusi watermark
Ldek=3;%Level dekomposisi Wavelet
jenis=2;%2:QIM-2
nbit=16;
nb=5;%bit kuantisasi QIM
mode=1;%0=qim,1=ss
alfa=0.1;

%%%%%%%%%%%%%%%%%%%
%Aktifkan yang sesuai Grup
%%%%%%%%%%%%%%%%%%%
%[1 (swt/dwt/lwt on) 1 (dst on) 1 (dct/fft on) 1/2 (svd/qr) 1 (cpt on)];
% sub_eksis=[2 0 0 1 1];%berhasil di semua subband
% sub_eksis=[2 1 0 1 1];%berhasil di semua subband (Grup 3)
% sub_eksis=[2 0 1 1 1];%berhasil di semua subband
% sub_eksis=[2 1 1 1 1];%berhasil di semua subband (Grup 1)
% sub_eksis=[2 0 0 2 1];%berhasil di subband 1
% sub_eksis=[2 1 0 2 1];%berhasil di semua subband (Grup 5)
 sub_eksis=[2 0 1 2 1];%berhasil di semua subband  (Grup 6)
% sub_eksis=[2 1 1 2 1];%berhasil di semua subband (Grup 2)

% alfass=0.00000000000001; %dwt alfass minimal ber =0, subbad 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%SILAKAN DIUBAH-UBAH%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subbad=1;%subband
alfass=0.001; 
B=16; %Resolusi audio BxB saat diubah ke 2D, atau 1 segmennya B^2 sampel
serangan=[8 3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


folderhost=[pwd '/host_audio/'];
folderwatermark=[pwd '/watermark/'];

%Membaca file audio host
[x,fs]=audioread([folderhost 'evangeline-matthew_sweet.wav']);
% [x,fs]=audioread(['output.wav']);

%Panjang sampel dalam 1 segmen
L=B^2;
% B=sqrt(L);
%Audio dibatasi sebanyak seluruh watermark yang disisipkan
x1=x(1:L*M*N,1);
% x1=0.1*ones(size(x1));
%Membaca file watermark
logo=imread([folderwatermark 'logo tel-u.png']);
%Mengubah ukuran watermark,
logor=imresize(logo(1:round(size(logo,1)*0.75),:,:),[M N]);
%Ubah ke grayscale
logogray=rgb2gray(logor);
%Ubah ke B/W
logobw=~im2bw(logogray,0.0000001);
%Ubah ke 1-D, panjang watermark=MN
logobw1d=double(logobw(:));
%Ubah watermark dari RZ ke NRZ
logonrz=2*logobw1d-1;

if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT
    %Generate pncode
    pn0=2*randi([0 1],B,B)-1;
    pn1=2*randi([0 1],B,B)-1;
elseif sub_eksis(1)==2 %DWT
    %Generate pncode
    pn0=2*randi([0 1],B/2,B/2)-1;
    pn1=2*randi([0 1],B/2,B/2)-1;
end

    if sub_eksis(4)==1 %SVD
% pn0=randi([0 1],B,B);
% pn1=randi([0 1],B,B);
[U0,S0,V0]=svd(pn0);
[U1,S1,V1]=svd(pn1);
    elseif sub_eksis(4)==2

        [U0,S0]=qr(pn0);
        [U1,S1]=qr(pn1);
    % elseif sub_eksis(4)==0
    %     S=Swdct2;
    end

% % Embedding dengan Metode SS
% for i=1:length(logonrz)
%     xw(L*i-L+1:L*i,1)=x1(L*i-L+1:L*i)+alfa*logonrz(i)*pn;
% end
% mean(mean((u*s*v.').*pn))
% mean(mean((u*s*v.').*pn1))
% mean(mean((u1*s1*v1.').*pn))
% mean(mean((u1*s1*v1.').*pn1))
% mean(mean((u*s*v.').*pn))
% mean(mean((u*s1*v.').*pn1))

%%%%%%%%%%%%%%EMBEDDING%%%%%%%%%%%%%%
for i=1:M*N
    %Segmentasi
    seg=alfa*x1(L*i-L+1:L*i,1);%segmentasi
%Wavelet Transform
    if sub_eksis(1)==1 %SWT
        Sw=swt(seg,Ldek,'haar');
    elseif sub_eksis(1)==0
        Sw=seg.';subbad=1;
    elseif sub_eksis(1)==2
        [A,D]=dwt(seg','haar');
        [LL,LH]=dwt(A,'haar');
        [HL,HH]=dwt(D,'haar');
        Sw=[LL;LH;HL;HH];
    end
%DST
    if sub_eksis(2)==0
        sdst=Sw(subbad,:);
    elseif sub_eksis(2)==1 %DST
        sdst=dst(Sw(subbad,:));
    end
    %DCT
    if sub_eksis(3)==0
        Swdct=sdst;
    elseif sub_eksis(3)==1 %DCT
        Swdct=dct(sdst);
    end
    %1D ke 2D
    if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT

        Swdct2=reshape(Swdct,[B B]);
    elseif sub_eksis(1)==2 %DWT
        Swdct2=reshape(Swdct,[B/2 B/2]);

    end
    %SVD atau QR
    if sub_eksis(4)==1 %SVD
        [U,S,V]=svd(Swdct2);
    elseif sub_eksis(4)==2
        [U,S]=qr(Swdct2);
    elseif sub_eksis(4)==0
        S=Swdct2;
    end
    if mode==0
    %untuk QIM
        %CPT dengan berbagai parameter
        if sub_eksis(5)==0 %Tanpa CPT
            swcpt(i)=S(1);
        elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
            %SVD-CPT
            [r, theta, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=r;
        elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=theta_deg;
        elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=theta_rad;
        elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=r;
        elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=rad2deg(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=rad2deg(phi);
        elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=(phi);
        end

        %QIM embedding
        % r_emb(i)=QIM_emb(r,logobw1d(i),nb,jenis);
        % theta_emb(i)=QIM_emb(theta_deg,logobw1d(i),nb,jenis);
        qim_emb(i)=QIM_emb(swcpt(i),logobw1d(i),nb,jenis);
        % qim_emb(i)=swcpt(i);
        wt(i,1)=QIM_ext(qim_emb(i),nb,jenis);
        Sdgab=S;
        if sub_eksis(5)==0 %Tanpa CPT
            Sdgab(1)=qim_emb(i);


        elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
            %SVD-CPT
            r=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r, theta);
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);

        elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r,theta);
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);
        elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r,rad2deg(theta));
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
            r=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, deg2rad(theta), phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
            phi=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, deg2rad(phi));
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
            phi=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        end

    elseif mode==1
        %Spread Spectrum
        if sub_eksis(4)==1 %SVD
            So(i).data={S};
            if logonrz(i)==1
                Sdgab=S+alfass*S1;
            else
                Sdgab=S+alfass*S0;
            end
            X00=(mean(mean(abs((U0*(Sdgab-cell2mat(So(i).data))/alfass*V0.').*pn0))));
            X11=(mean(mean(abs((U1*(Sdgab-cell2mat(So(i).data))/alfass*V1.').*pn1))));
            if X00>=X11
                wt(i,1)=0;
            else
                wt(i,1)=1;
            end
        elseif sub_eksis(4)==2 %QR
            So(i).data={S};
            if logonrz(i)==1
                Sdgab=S+alfass*S1;
            else
                Sdgab=S+alfass*S0;
            end
            X00=(mean(mean(abs((U0*(Sdgab-cell2mat(So(i).data))/alfass).*pn0))));
            X11=(mean(mean(abs((U1*(Sdgab-cell2mat(So(i).data))/alfass).*pn1))));
            if X00>=X11
                wt(i,1)=0;
            else
                wt(i,1)=1;
            end

        elseif sub_eksis(4)==0 %

        end

    end


    % %QIM extraction
    % wt(i,1)=QIM_ext(Semb(i),nb,jenis);
    %Penggabungan nilai Singular yg disisipkan dg yg tdk
    if sub_eksis(4)==1 %ISVD
        Sisvd=U*Sdgab*V.';
    elseif sub_eksis(4)==2%QR Rec.
        Sisvd=U*Sdgab;
    elseif sub_eksis(4)==0
        Sisvd=Sdgab;

    end
    %2D ke 1D
    if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT

        Sisvd1=reshape(Sisvd,[(B)^2 1]);
    elseif sub_eksis(1)==2  %DWT
        Sisvd1=reshape(Sisvd,[(B/2)^2 1]);
    end
    %IDCT
    %Menggabungkan subband yg disisipkan dg yg tdk disisipkan
    Swgab=Sw;
    if sub_eksis(3)==0
        s_idst=Sisvd1;
    elseif sub_eksis(3)==1
        s_idst=idct(Sisvd1);
    end
    
    %IDST
    if sub_eksis(2)==0
        %IDST
        Swgab(subbad,:)=s_idst;

    elseif sub_eksis(2)==1
        %IDST
        Swgab(subbad,:)=idst(s_idst);
    end

    %Wavelet Reconstruction
    if sub_eksis(1)==1
        %ISWT
        st=iswt(Swgab,'haar');
    elseif sub_eksis(1)==0
        st=Swgab;
    elseif sub_eksis(1)==2
        cA=idwt(Swgab(1,:),Swgab(2,:),'haar');
        cD=idwt(Swgab(3,:),Swgab(4,:),'haar');
        st=idwt(cA,cD,'haar');
    end

    %Penggabungan Segmen
    xw(L*i-L+1:L*i,1)=st/alfa;
end

hasil.snr=10*log10((mean(x1.^2))/(mean((x1-xw).^2)));
hasil.odg=hitungodgsegmen(1,xw,x1,fs,nbit);

% xw=0.9*xw;

hasil.payload=fs/L;
audiowrite('temphost.wav',xw,fs,'BitsPerSample',nbit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Serangan Audio Watermarking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% serangan=[0 0];
% xwn=audio_attack(xw,fs,serangan,snrnoise,fr,A,f);
allattack_audio_stirmark(serangan,'temphost.wav','temp_attack.wav',nbit)
% xwn=xw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xwn,fs]=audioread(['temp_attack.wav']);

%%%%%%%%%%%%%EXTRACTION%%%%%%%%%%%%%%
for i=1:M*N
    % for i=1:1

    %Segmentasi
    seg=alfa*xwn(L*i-L+1:L*i,1);
    %SWT
    if sub_eksis(1)==1 %SWT
        Sw=swt(seg,Ldek,'haar');
    elseif sub_eksis(1)==0
        Sw=seg.';
    elseif sub_eksis(1)==2
        [A,D]=dwt(seg','haar');
        [LL,LH]=dwt(A,'haar');
        [HL,HH]=dwt(D,'haar');
        Sw=[LL;LH;HL;HH];

    end
    %DST
    if sub_eksis(2)==0
        sdst=Sw(subbad,:);
    elseif sub_eksis(2)==1 %DST
        sdst=dst(Sw(subbad,:));
    end

    %DCT
    if sub_eksis(3)==0
        Swdct=sdst;
    elseif sub_eksis(3)==1 %DCT
        Swdct=dct(sdst);
    end
    %1D ke 2D
    if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT
        Swdct2=reshape(Swdct,[B B]);
    elseif sub_eksis(1)==2  %DWT
        Swdct2=reshape(Swdct,[B/2 B/2]);
    end

    %SVD atau QR
    if sub_eksis(4)==1 %SVD
        [U,S,V]=svd(Swdct2);
    elseif sub_eksis(4)==2
        [U,S]=qr(Swdct2);
    elseif sub_eksis(4)==0
        S=Swdct2;

    end
    if mode==0 %QIM
        %CPT dengan berbagai parameter
        if sub_eksis(5)==0 %Tanpa CPT
            swcptr(i)=S(1);
        elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
            %SVD-CPT
            [r, theta, theta_deg] = ctp(S(1),S(2));
            swcptr(i)=r;
        elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcptr(i)=theta_deg;
        elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcptr(i)=theta_rad;
        elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcptr(i)=r;
        elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcptr(i)=rad2deg(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcptr(i)=(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcptr(i)=rad2deg(phi);
        elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcptr(i)=(phi);
        end
        %QIM extraction
        % wt(i,1)=QIM_ext(rext(i),nb,jenis);
        % wt(i,1)=QIM_ext(theta_deg(i),nb,jenis);
        wtt(i,1)=QIM_ext(swcptr(i),nb,jenis);
    elseif mode==1 %Spread Spectrum
        if sub_eksis(4)==1 %SVD
            X0=(mean(mean(abs((U0*(S-cell2mat(So(i).data))/alfass*V0.').*pn0))));
            X1=(mean(mean(abs((U1*(S-cell2mat(So(i).data))/alfass*V1.').*pn1))));
            if X0>=X1
                wtt(i,1)=0;
            else
                wtt(i,1)=1;
            end
        elseif sub_eksis(4)==2 %QR
            X0=(mean(mean(abs((U0*(S-cell2mat(So(i).data))/alfass).*pn0))));
            X1=(mean(mean(abs((U1*(S-cell2mat(So(i).data))/alfass).*pn1))));
            if X0>=X1
                wtt(i,1)=0;
            else
                wtt(i,1)=1;
            end
        elseif sub_eksis(4)==0 %

        end


    end



% figure(1),clf
% plot(diag(cell2mat(So(i).data)),'b-'),hold on,grid on
% plot(diag(Sdgab),'g-')
% plot(diag(S),'r-'),
% legend('Asli','Sw','Swr')

% figure(1),clf
% plot(qim_emb,'b-'),hold on,grid on
% plot(swcpt,'r-','linewidth',2)
% plot(swcptr,'g-')
% legend('out qim','before qim','in rx')
% 
% continue
end

hasil.ber0=mean(mean(abs(wt-double(logobw1d))));
hasil.ber1=mean(mean(abs(wtt-double(logobw1d))));
hasil

logo2d=reshape(logobw1d,[M N]);
logo2dt=reshape(wtt,[M N]);
figure(1),clf
subplot(121),imshow(uint8(255*logo2d)),title('(a)')
subplot(122),imshow(uint8(255*logo2dt)),title('(b)')
% figure(1),clf
% plot(diag(cell2mat(So(1).data)),'b-'),hold on,grid on
% plot(diag(Sdgab),'g-','LineWidth',2)
% plot(diag(S),'r-'),
% legend('Asli','Sw','Swr')

% return
% figure(1),clf
% plot(qim_emb,'b-'),hold on,grid on
% plot(swcpt,'r-','linewidth',2)
% plot(swcptr,'g-')
% legend('out qim','before qim','in rx')


return


