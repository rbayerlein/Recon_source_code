function [lm_all, a_all, m_all] = combine_all_8beds(data_path)

lm_all=[];
a_all=[];
m_all=[];

for n=1:8, n

%fname=sprintf('04-14-2014/Grp0_Frm0_Sub%d_pp_tofMask_tof_corr.prl',n);
%cmd=sprintf('../lmdata_toshiba recon150.cfg %s', fname);
%unix(cmd);

%lm=touch2('my.lm','int16');
lm = touch2(sprintf('%s/bed_%d.lm', data_path, n), 'int16');

lm=reshape(lm, 5, length(lm)/5);

%flip now
lm(2,:) = 48 - lm(2,:) - 1;
lm(4,:) = 48 - lm(4,:) - 1;

% shift
lm(2,:) = lm(2,:) + 24 * (n-1);
lm(4,:) = lm(4,:) + 24 * (n-1);

a = touch(sprintf('%s/bed_%d.add_fac', data_path, n));
m = touch(sprintf('%s/bed_%d.mul_fac', data_path, n));
%a=touch('my.add_fac');
%m=touch('my.mul_fac');

lm_all = [lm_all, lm];
a_all=[a_all; a];
m_all=[m_all; m];

end
