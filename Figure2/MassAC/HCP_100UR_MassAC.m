clear

addpath(genpath('/home/wmrnaq/bin/xDF/'))
addpath(genpath('/home/wmrnaq/bin/HetBiv/'))

DirID = 'LR';

load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

i = str2double(getenv('SGE_TASK_ID'));

SubID = SubList{i};

WhereFrom = ['/storage/essicd/data/HCP/Soroosh/HCP_100Unrelated/' SubID '_' DirID '/'];

[Y,Y_Stat] = CleanNIFTI([WhereFrom 'rfMRI_REST1_' DirID '_hp2000_clean.nii.gz'],'demean');
T = size(Y,2);
I = size(Y,1);

%Because we don't have enough memory on the damning damn Buster!
Smpl = 10000;
for i = 1:round(I./Smpl)
        disp(['Y:' num2str(size(Y))])
        idx_tmp = (i-1)*Smpl+1:i*Smpl;

	if max(idx_tmp)>I
		idx_tmp = (i-1)*Smpl+1:I;
	end

        AC(idx_tmp,:) = AC_fft(Y(1:Smpl,:),T);
        Y(1:Smpl,:)=[];
        
        if size(Y,1)<Smpl
            disp(['There was ' num2str(size(Y,1)) ' remainders.'])
	    AC = [AC; AC_fft(Y,T)];
	    break
        end
        
        disp(['AC:' num2str(size(AC))])
end
size(AC), size(Y)
%---------------------------------------------------------------


idx = 1:Y_Stat.OrigDim(1);
idx(Y_Stat.Removables) = [];
YImg = zeros(Y_Stat.OrigDim);
YImg(idx,:) = AC;
YImg = reshape(YImg,[Y_Stat.ImgDim]); %Okay! this is clean handsome 3D volume!

ImgObj = Y_Stat.Obj;
ImgObj.img = YImg;

disp(['Saving the AC image...'])

Where2 = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_MassAC/' SubID '_' DirID '/']

if exist(Where2,'dir')~=7; mkdir(Where2); end

save_untouch_nii(ImgObj,[Where2 'AC_' SubID '_' DirID '.nii.gz']); %this is the cleanest version, ready to be parcelled. 
disp('Done!')
