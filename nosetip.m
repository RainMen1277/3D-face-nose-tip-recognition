clc;
read_file = importdata('C:\myfiles\xyz\test1\test3\3.xyz'); %读取文件
file_x = read_file(:,1);
file_y = read_file(:,2);
file_z = read_file(:,3);
set(figure,'name','点云');
plot3(file_x,file_y,file_z,'.')
[X,Y,Z] = griddata(file_x,file_y,file_z,linspace(min(file_x),max(file_x),200)',linspace(min(file_y),max(file_y),20),'cubic'); %插值
set(figure,'name','网格');
surf(Z);  %显示三维曲面
nose = []; %鼻尖候选点
x_gaos = []; %鼻尖置信点
II = [];
% 求Y轴的切片有效范围
for j_new = 1:20 %鲁莽的直接选择切片高度
    xz = Z(j_new,:); %获取水平切片
    %画圆
    I = find(~isnan(xz)); %获取水平切片非空位置
    [x_s,~] = min(I); %水平切片有效起始位置
    [x_e,~] = max(I); %水平切片有效终止位置
    %r = round((x_e-x_s)/5); %自适应因切圆半径(效果并不好)
    r = 30; %指定切圆的半斤，效果更好
    h_gaos = []; %切片中每个点作为圆心和切片形成的交点所形成的三角形的高集合
    for x1 = x_s+r:x_e-r %圆心的可取值范围
        if x1>1 %避免水平切片太小
            z1 = xz(1,x1); %切片中对应z的坐标
            i = x1-r:x1+r; %圆心的取值范围
            y2 = sqrt(r.^2-(i-x1).^2)+z1; %切圆的z坐标（上半圆）
            y = xz(1,i);
            cy = y2-y;
            pos = cy>0;
            neg = cy<=0;
            %确定变号位置
            fro = diff(pos)~=0; %变号的前导位置
            rel = diff(neg)~=0; %变号的尾巴位置
            zpf = find(fro==1); %记录索引
            zpr = find(rel==1)+1; %记录索引
            zpfr = [zpf,zpr];
            x0 = (i(zpr).*(y2(zpf)-y(zpf))-i(zpf).*(y2(zpr)-y(zpr)))./(y(zpr)+y2(zpf)-y(zpf)-y2(zpr));
            y0 = y(zpf)+(x0-i(zpf)).*(y(zpr)-y(zpf))./(i(zpr)-i(zpr)-i(zpf));
            x0 = [x0 x0].';
            y0 = [y0 y0].';
            jie = unique([x0,y0],'rows'); %得到切片与下半圆的交点集合
            [mm,nn] = size(jie);
            %fprintf('解的个数：');
            %disp(mm);
            if mm == 2;
                jie1 = jie(1,:);
                jie2 = jie(2,:);
                if jie1(:,1)<x1 && jie2(:,1)>x1 %确保两个交点在圆心的两侧
                    l1 = ((jie1(:,1)-x1)^2+(jie1(:,2)-z1)^2)^(1/2);
                    l2 = ((jie2(:,1)-x1)^2+(jie2(:,2)-z1)^2)^(1/2);
                    l3 = ((jie1(:,1)-jie2(:,1))^2+(jie1(:,2)-jie2(:,2))^2)^(1/2); %内接三角形的底为
                    p = (l1+l2+l3)/2;
                    s = sqrt(p*(p-l1)*(p-l2)*(p-l3)); %内接三角形的面积
                    h_gao = 2*s/l3;  %内接三角形的高为
                    h_gaos = [h_gaos;[x1,h_gao]]; %形成的高的集合
                end
            end
        end
    end
    if length(h_gaos)>0
        [m,p] = max(h_gaos(:,2));
        x_gao = h_gaos(p,1); %得到最大的高的圆心x坐标
        x_gaos = [x_gaos;[j_new,x_gao,m]]; %切片Y坐标，圆心X坐标，最大三角形的高
    end
end
[max_h_1,max_p_1] = max(x_gaos(:,3));
j_max_1 = x_gaos(max_p_1,1);
j_2_s = Y((j_max_1-1),1);
j_2_e = Y((j_max_1+1),1);

number = size(read_file,1);
new_file = [];
for i = 1:number
    if j_2_s < read_file(i,2) && read_file(i,2) < j_2_e
        new_file = [new_file;read_file(i,:)];
    end
end
new_x = new_file(:,1);
new_y = new_file(:,2);
new_z = new_file(:,3);
[Xx,Yy,Zz] = griddata(new_x,new_y,new_z,linspace(min(new_x),max(new_x),200)',linspace(min(new_y),max(new_y),40),'cubic'); %插值

 for j_new = 1:40 %鲁莽的直接选择切片高度
    xz = Zz(j_new,:); %获取水平切片
    %画圆
    I = find(~isnan(xz)); %获取水平切片非空位置
    [x_s,~] = min(I); %水平切片有效起始位置
    [x_e,~] = max(I); %水平切片有效终止位置
    %r = round((x_e-x_s)/5); %自适应因切圆半径(效果并不好)
    r = 30; %指定切圆的半斤，效果更好
    h_gaos = []; %切片中每个点作为圆心和切片形成的交点所形成的三角形的高集合
    for x1 = x_s+r:x_e-r %圆心的可取值范围
        if x1>1 %避免水平切片太小
            z1 = xz(1,x1); %切片中对应z的坐标
            i = x1-r:x1+r; %圆心的取值范围
            y2 = sqrt(r.^2-(i-x1).^2)+z1; %切圆的z坐标（上半圆）
            y = xz(1,i);
            cy = y2-y;
            pos = cy>0;
            neg = cy<=0;
            %确定变号位置
            fro = diff(pos)~=0; %变号的前导位置
            rel = diff(neg)~=0; %变号的尾巴位置
            zpf = find(fro==1); %记录索引
            zpr = find(rel==1)+1; %记录索引
            zpfr = [zpf,zpr];
            x0 = (i(zpr).*(y2(zpf)-y(zpf))-i(zpf).*(y2(zpr)-y(zpr)))./(y(zpr)+y2(zpf)-y(zpf)-y2(zpr));
            y0 = y(zpf)+(x0-i(zpf)).*(y(zpr)-y(zpf))./(i(zpr)-i(zpr)-i(zpf));
            x0 = [x0 x0].';
            y0 = [y0 y0].';
            jie = unique([x0,y0],'rows'); %得到切片与下半圆的交点集合
            [mm,nn] = size(jie);
            %fprintf('解的个数：');
            %disp(mm);
            if mm == 2;
                jie1 = jie(1,:);
                jie2 = jie(2,:);
                if jie1(:,1)<x1 && jie2(:,1)>x1 %确保两个交点在圆心的两侧
                    l1 = ((jie1(:,1)-x1)^2+(jie1(:,2)-z1)^2)^(1/2);
                    l2 = ((jie2(:,1)-x1)^2+(jie2(:,2)-z1)^2)^(1/2);
                    l3 = ((jie1(:,1)-jie2(:,1))^2+(jie1(:,2)-jie2(:,2))^2)^(1/2); %内接三角形的底为
                    p = (l1+l2+l3)/2;
                    s = sqrt(p*(p-l1)*(p-l2)*(p-l3)); %内接三角形的面积
                    h_gao = 2*s/l3;  %内接三角形的高为
                    h_gaos = [h_gaos;[x1,h_gao]]; %形成的高的集合
                end
            end
        end
    end
    if length(h_gaos)>0
        [m,p] = max(h_gaos(:,2));
        x_gao = h_gaos(p,1); %得到最大的高的圆心x坐标
        x_gaos = [x_gaos;[j_new,x_gao,m]]; %切片Y坐标，圆心X坐标，最大三角形的高
    end
end
ans = sortrows(x_gaos,-3); %按照高的大小顺序排列（倒序）
nose_cut = ans(1:20,:); %取前20个最大的高对应的圆心坐标
for mq = 1:20
    nose = [nose;[Xx(1,nose_cut(mq,2)),Yy(nose_cut(mq,1),1),Zz(nose_cut(mq,1),nose_cut(mq,2))]];
end
% fprintf('鼻尖置信点的维度：');
% disp(x_gaos);
% fprintf('列出置信点：');
% disp(nose_cut);
% fprintf('在原坐标系中显示出鼻尖置信点的位置：');
% disp(nose);
hold on;
nose_x = nose(:,1);
nose_y = nose(:,2);
nose_z = nose(:,3);
plot3(nose_x,nose_y,nose_z,'o');

% RANSAC直线拟合
data=[nose(:,1)';nose(:,2)'];
iter = 500; 
number = size(data,2); % 总点数
bestParameter1=0; bestParameter2=0; % 最佳匹配的参数
sigma = 1;
pretotal=0;     %符合拟合模型的数据的个数
for i=1:iter
    idx = randperm(number,2);  %随机选择两个点
    sample = data(:,idx); 
    line = zeros(1,3);%拟合直线方程 y=kx+b
    x = sample(1, :);
    y = sample(2, :);
    if x(1)==x(2)
        a = 1;
        b = 0;
        c = -x(1);
        line = [a b c];
    else
        k=(y(1)-y(2))/(x(1)-x(2)); %直线斜率
        a = sqrt(1-1/(1+k^2));
        b = sqrt(1-a^2);
        if k > 0
            b = -b;
        end
        c = -a*x(1)-b*y(1);
        line = [a b c];
    end
    mask=abs(line*[data; ones(1,size(data,2))]);    %求每个数据到拟合直线的距离
    total=sum(mask<sigma);              %计算数据距离直线小于一定阈值的数据的个数
    if total>pretotal            %找到符合拟合直线数据最多的拟合直线
        pretotal=total;
        bestline=line;          %找到最好的拟合直线
    end  
end
%显示符合最佳拟合的数据
mask=abs(bestline*[data; ones(1,size(data,2))])<sigma;    
hold on;
k=1;
nose_z = 1000;
for i=1:length(mask)
    if mask(i)
        if nose(i,3) < nose_z
            nose_z = nose(i,3);
            nose_i = i;
        end
    end
end
hold on;
nose_max_x = nose(nose_i,1);
nose_max_y = nose(nose_i,2);
nose_max_z = nose(nose_i,3);
plot3(nose_max_x,nose_max_y,nose_max_z,'g*');