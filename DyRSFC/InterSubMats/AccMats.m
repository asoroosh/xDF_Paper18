clear

nRlz = 5000;

for i = 1:nRlz

	if ~mod(i,100); disp(num2str(i)); end;
	
        R = load(['/storage/essicd/data/HCP/Soroosh/DyConn/InterSubMats/YEO/InterSubMats_YEO_' num2str(i) '.mat']);
	
        mat_naive(:,:,i) = R.Mat.mat_naive;
	mat_xdff(:,:,i)  = R.Mat.mat_xdff; 
	mat_xdf(:,:,i)   = R.Mat.mat_xdf;
	mat_wn(:,:,i)    = R.Mat.mat_wn;
	mat_pw(:,:,i)    = R.Mat.mat_pw;

end


save(['R/HCP_InterSub_DyConn_YEO.mat'],'mat_naive','mat_xdff','mat_xdf','mat_wn','mat_pw')
