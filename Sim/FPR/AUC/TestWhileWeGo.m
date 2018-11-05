clear 

for T = [100 200 600 1200]
	for i=1:2000
		WhereFrom = ['R/Sen_Spec_t' num2str(T) '_118528_' num2str(i) '_TVOff_r9-2.mat'];
%		WhereFrom
		if ~exist(WhereFrom,'file');
			disp(num2str(i))
			continue;
		end

		R=load(WhereFrom);
		Acc(i,:) = R.Acc;
		Sen(i,:) = R.Sen;
		Spc(i,:) = R.Spc;
	end
	mean(Acc)
	mean(Sen)
	mean(Spc)

clear Acc Sen Spc
end
