gunzip('T1_r.nii.gz'); %% T1_r.nii.gz --> T1_r.nii

start_seg('T1_r.nii',[]); %% segmentation T1_r.nii가 있는 폴더에서 실행

segTouchup('T1_r.nii',[]); %% segmentation 된 이미지 보정

electrode = '10_20'; %%10_10 =81, 아니면 18

if electrode == '10_10'
    load('81options.mat')
    hdrinfo=electrodePlacement('T1_r.nii','T1_r.nii',[],{'F3';'F4'},options,int2str(1));
    num_el=81;
    reference=-3;
else
    load('18options.mat')
    hdrinfo=electrodePlacement('T1_r.nii','T1_r.nii',[],{'F3';'F4'},options,int2str(1));
    num_el=18;
    reference=-3;
end

[node,elem,face] = meshByIso2mesh('T1_r.nii','T1_r.nii',[],options.meshOpt,hdrinfo,int2str(1)); %% mesh 생성

%% solve하는 과정
prepareForGetDP('T1_r.nii',node,elem,{'F3';'F4'},int2str(1),num_el); %% solve하기 전에 하는거

if isequal(app.electrode, '10_10')
    indElecSolve = 1:num_el;
    temp=zeros(num_el,1);
    if ~exist('T1_r_81_e.pos')
        for i=1:81
        temp=zeros(81,1);
        temp(i)=1;temp(81+reference)=-1;
        if i==81+reference
        else
            solveByGetDP('T1_r.nii',temp,options.conductivities,indElecSolve,int2str(i),'');
        end
        end
    end
else
    indElecSolve = 1:num_el;
    temp=zeros(num_el,1);
     if ~exist('T1_r_77_e.pos')
        for i=[1,3,14,16,18,20,22,36,38,40,42,44,56,58,60,62,64,77]
        temp=zeros(81,1);
        temp(i)=1;temp(81+reference)=-1;
        if i==81+reference
        else
            solveByGetDP('T1_r.nii',temp,options.conductivities,indElecSolve,int2str(i),'');
        end
        end
     end
end
%% optimization matrix 만드는거
if isequal(electrode, '10_10')
    E_all={};
    tt=0;
    for i=1:81
        if i==81+reference
        else
            tt=tt+1;
            temp_name=['T1_r_' int2str(i) '_e.pos'];
            E_all{tt}=dlmread(temp_name);
        end
    end
    load('opt_matrix81.mat')
else
    E_all={};
    tt=0;
    for i=[1,3,14,16,18,20,22,36,38,40,42,44,56,58,60,62,64,77]
        if i==81+reference
        else
            tt=tt+1;
            temp_name=['T1_r_' int2str(i) '_e.pos'];
            E_all{tt}=dlmread(temp_name);
        end
    end
    hey = load('opt_matrix18.mat').bb;
end

%% 특정 layer에서 좌표 받아와 target 영역으로 지정하고 해당 영역의 TI effect를 최적화화는 알고리즘

layer=1; %% 
brain=unique(elem(elem(:,5)==layer,1:4));
target=[137.1 141 148.5]; %% target1=52.49,103.5,168; target2=52.5 67.24 154.8; target3=94.12 114.5 134.2 
%%51.37 152.1 142.1

nl=find(abs(node(brain,1)-target(1))<3&abs(node(brain,2)-target(2))<3&abs(node(brain,3)-target(3))<3);
tic
E_temp=zeros(length(bb),9);
E_dot=zeros(length(bb),1);
for i=1:length(bb)
    
    E1=E_all{1,bb(i,1)}(2:end,2:end)-E_all{1,bb(i,2)}(2:end,2:end);
    E2=E_all{1,bb(i,3)}(2:end,2:end)-E_all{1,bb(i,4)}(2:end,2:end);

    E=[];

    E_dot(i,k)=sum(E1(brain(nl),1).*E2(brain(nl),1)+E1(brain(nl),2).*E2(brain(nl),2)...
        +E1(brain(nl),3).*E2(brain(nl),3));
end
b1=find(E_dot==max(E_dot,[],'all'));

for k=1:9
        
        E1=E_all{1,bb(b1,1)}(2:end,2:end)-E_all{1,bb(b1,2)}(2:end,2:end);
        E2=E_all{1,bb(b1,3)}(2:end,2:end)-E_all{1,bb(b1,4)}(2:end,2:end);

        E=[];

        for j=1:length(nl)
           E(j)=calc_TI(0.2*k*E1(brain(nl(j)),:),(2-0.2*k)*E2(brain(nl(j)),:));
        end
        E_temp(k)=sum(E);
end

toc
[a1]=find(E_temp==max(E_temp,[],'all'));

% optldx=3;

optldx=b1;
E1=E_all{1,bb(optldx,1)}(2:end,2:end)-E_all{1,bb(optldx,2)}(2:end,2:end);
E2=E_all{1,bb(optldx,3)}(2:end,2:end)-E_all{1,bb(optldx,4)}(2:end,2:end);
E=[];
for j=1:length(node)
  E(j)=calc_TI(0.2*a1*E1(j),(2-0.2*a1)*E2(j));
end


% scatter3(node(brain,1),node(brain,2),node(brain,3),5,E(brain)')
% hold on
% plot3(node(brain(nl),1),node(brain(nl),2),node(brain(nl),3),'o','Color','r','MarkerSize',10,'MarkerFaceColor',[1 0 0])
% caxis([0 0.3]);

%% TI 최적화 결과
plotsurf([node(:,1:3),E'],face(face(:,4)==1,1:4),'EdgeColor','None')
caxis([0 0.3]);

light("Style","local","Position",[200 0 0]);
light("Style","local","Position",[-200 0 0]);
light("Style","local","Position",[0 200 0]);
light("Style","infinite","Position",[0 0 200]);
view([-1 1 1])
axis off
set(gca, 'Position', [-0.25, -0.2, 1.5, 1.5]);
set(gcf,'color','w');
