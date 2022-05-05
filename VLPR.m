clear ;
close all;

Scolor = imread('Original image.jpg');

Sgray = rgb2gray(Scolor);
figure,imshow(Scolor),title('Original color image');
figure,imshow(Sgray),title('Original black-and-white image');

s=strel('disk',13);
Bgray=imopen(Sgray,s);
figure,imshow(Bgray);title('Background image');

Egray=imsubtract(Sgray,Bgray);
figure,imshow(Egray);title('Enhance black-and-white image');

fmax1=double(max(max(Egray)));
fmin1=double(min(min(Egray)));
level=(fmax1-(fmax1-fmin1)/3)/255;
bw22=im2bw(Egray,level);
bw2=double(bw22);

figure,imshow(bw2);title('Image binarization');
grd=edge(bw2,'canny')
figure,imshow(grd);title('Image edge extraction');
bg1=imclose(grd,strel('rectangle',[5,19]));
figure,imshow(bg1);title('Image closed operation[5,19]');
bg3=imopen(bg1,strel('rectangle',[5,19]));
figure,imshow(bg3);title('Image open operation[5,19]');
bg2=imopen(bg3,strel('rectangle',[19,1]));
figure,imshow(bg2);title('Image open operation[19,1]');

[L,num] = bwlabel(bg2,8);
Feastats = regionprops(L,'basic');
Area=[Feastats.Area];
BoundingBox=[Feastats.BoundingBox];
RGB = label2rgb(L, 'spring', 'k', 'shuffle'); 
figure,imshow(RGB);title('Color mark of image');
lx=0;
for l=1:num
    width=BoundingBox((l-1)*4+3);
    hight=BoundingBox((l-1)*4+4);
    if (width>98 & width<160 & hight>25 & hight<50)
        lx=lx+1;
        Getok(lx)=l;
    end
end
for k= 1:lx
    l=Getok(k);    
    startcol=BoundingBox((l-1)*4+1)-2;
    startrow=BoundingBox((l-1)*4+2)-2;
    width=BoundingBox((l-1)*4+3)+8;
    hight=BoundingBox((l-1)*4+4)+2;
    rato=width/hight;
    if rato>2 & rato<4   
        break;
    end
end
sbw1=bw2(startrow:startrow+hight,startcol:startcol+width-1); 
subcol1=Sgray(startrow:startrow+hight,startcol:startcol+width-1);
figure,subplot(2,1,1),imshow(subcol1);title('Gray subgraph of license plate');
subplot(2,1,2),imshow(sbw1);title('Binary subgraph of license plate');

histcol1=sum(sbw1);      
histrow=sum(sbw1');      
figure,subplot(2,1,1),bar(histcol1);title('Vertical projection (including border)');
subplot(2,1,2),bar(histrow);     title('Horizontal projection (including border)');
figure,subplot(2,1,1),bar(histrow);     title('Horizontal projection (including border)');
subplot(2,1,2),imshow(sbw1);title('Binary subgraph of license plate');

meanrow=mean(histrow);
minrow=min(histrow);
levelrow=(meanrow+minrow)/2;
count1=0;
l=1;
for k=1:hight
    if histrow(k)<=levelrow                             
        count1=count1+1;                                
    else 
        if count1>=1
            markrow(l)=k;
            markrow1(l)=count1;
            l=l+1;
        end
        count1=0;
    end
end
markrow2=diff(markrow);
[m1,n1]=size(markrow2);
n1=n1+1;
markrow(l)=hight;
markrow1(l)=count1;
markrow2(n1)=markrow(l)-markrow(l-1);
l=0;
for k=1:n1
    markrow3(k)=markrow(k+1)-markrow1(k+1);
    markrow4(k)=markrow3(k)-markrow(k);
    markrow5(k)=markrow3(k)-double(uint16(markrow4(k)/2));
end 


[m2,n2]=size(sbw1);
[m1,n1]=size(markrow4);
maxw=max(markrow4);
if markrow4(1) ~= maxw
    ysite=1;
    k1=1;
    for l=1:n2
    for k=1:markrow3(ysite)
        if sbw1(k,l)==1
            xdata(k1)=l;
            ydata(k1)=k;
            k1=k1+1;
            break;
        end
    end
    end
else  
    ysite=n1;
    if markrow4(n1) ==0
        if markrow4(n1-1) ==maxw
           ysite= 0; 
       else
           ysite= n1-1;
       end
    end
    if ysite ~=0
        k1=1;
        for l=1:n2
            k=m2;
        while k>=markrow(ysite) 
                if sbw1(k,l)==1
                    xdata(k1)=l;
                    ydata(k1)=k;
                    k1=k1+1;
                    break;
                end
                k=k-1;
            end
        end
    end
end       

fresult = fit(xdata',ydata','poly1');   
p1=fresult.p1;
angle=atan(fresult.p1)*180/pi; 
subcol = imrotate(subcol1,angle,'bilinear','crop'); 
sbw = imrotate(sbw1,angle,'bilinear','crop');
figure,subplot(2,1,1),imshow(subcol);title('Gray subgraph of license plate');
subplot(2,1,2),imshow(sbw);title('');
title(['Rotation angle of license plate: ',num2str(angle),'¶È'] ,'Color','r');

histcol1=sum(sbw); 
histrow=sum(sbw'); 
figure,subplot(2,1,1),bar(histcol1);title('Vertical projection (after rotation)');
subplot(2,1,2),bar(histrow);     title('Horizontal projection (after rotation)');
figure,subplot(2,1,1),bar(histrow);     title('Horizontal projection (after rotation)');
subplot(2,1,2),imshow(sbw);title('Binary subgraph of license plate (after rotation)');

maxhight=max(markrow2);
findc=find(markrow2==maxhight);
rowtop=markrow(findc);
rowbot=markrow(findc+1)-markrow1(findc+1);
sbw2=sbw(rowtop:rowbot,:);  
maxhight=rowbot-rowtop+1;   

histcol=sum(sbw2);  
figure,subplot(2,1,1),bar(histcol);title('Vertical projection (after removing the horizontal border)');%
subplot(2,1,2),imshow(sbw2); 
title(['Character height of license plate£º ',int2str(maxhight)],'Color','r');

meancol=mean(histcol);
mincol=min(histcol);
levelcol=(meancol+mincol)/4;
count1=0;
l=1;
for k=1:width
    if histcol(k)<=levelcol 
        count1=count1+1;
    else 
        if count1>=1
            markcol(l)=k; 
            markcol1(l)=count1; 
            l=l+1;
        end
        count1=0;
    end
end
markcol2=diff(markcol);
[m1,n1]=size(markcol2);
n1=n1+1;
markcol(l)=width;
markcol1(l)=count1;
markcol2(n1)=markcol(l)-markcol(l-1);

l=0;
for k=1:n1
    markcol3(k)=markcol(k+1)-markcol1(k+1);
    markcol4(k)=markcol3(k)-markcol(k); 
    markcol5(k)=markcol3(k)-double(uint16(markcol4(k)/2));
end 
markcol6=diff(markcol5); 
maxs=max(markcol6); 
findmax=find(markcol6==maxs);
markcol6(findmax)=0;
maxwidth=max(markcol6);

l=1;
[m2,n2]=size(subcol);
figure;
for k=findmax-1:findmax+5
        cleft=markcol5(k)-maxwidth/2;
        cright=markcol5(k)+maxwidth/2-2;
        if cleft<1
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2
            cright=n2;
            cleft=n2-maxwidth;
        end
        SegGray=sbw(rowtop:rowbot,cleft:cright);
        SegBw1=sbw(rowtop:rowbot,cleft:cright);
        SegBw2 = imresize(SegBw1,[22 14]);
        subplot(2,n1,l),imshow(SegGray);
        if l==7
            title(['Character width of license plate£º ',int2str(maxwidth)],'Color','r');
        end
        subplot(2,n1,n1+l),imshow(SegBw2);               
        fname=strcat('E:\Matlab\bin',int2str(k),'.jpg');
        imwrite(SegBw2,fname,'jpg') 
        l=l+1;
end

liccode=char(['0':'9' 'A':'Z' '´¨¶õ¸Ó¹ð»¦¼ª½ò½ú¾©ÁÉÂ³ÃÉÃöÄþÉÂËÕÍîÏæÔ¥ÔÁÕã']); 
SubBw2=zeros(22,14);
l=1;
[m2,n2]=size(sbw);
for k=findmax-1:findmax+5
       cleft=markcol5(k)-maxwidth/2;
        cright=markcol5(k)+maxwidth/2-2;
        if cleft<1
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2
            cright=n2;
            cleft=n2-maxwidth;  end
        SegBw1=sbw(rowtop:rowbot,cleft:cright);
        SegBw2 = imresize(SegBw1,[22 14]);
        if l==1              
            kmin=37;
            kmax=45;
        elseif l==2             
            kmin=11;
            kmax=36;
        elseif l>=3 & l<=5    
            kmin=1;
            kmax=36;
        else                    
            kmin=1;
            kmax=10;
        end
        for k2=kmin:kmax
            fname=strcat('E:\Matlab\bin\char\',liccode(k2),'.bmp');
            SamBw2 = imread(fname);           
            for  i=1:22
                for j=1:14
                    SubBw2(i,j)=SegBw2(i,j)-SamBw2(i,j);
                end
            end
            Dmax=0;
            for k1=1:22
                for l1=1:14
                    if ( SubBw2(k1,l1) > 0 | SubBw2(k1,l1) <0 )
                        Dmax=Dmax+1;
                    end
                end
            end
            Error(k2)=Dmax;
        end
        Error1=Error(kmin:kmax);
        MinError=min(Error1);
        findc=find(Error1==MinError);
        RegCode(l*2-1)=liccode(findc(1)+kmin-1);
        RegCode(l*2)=' ';
        l=l+1;
end
title (['Identify the license plate number:', RegCode],'Color','r');