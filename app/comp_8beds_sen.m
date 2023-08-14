function s0 = comp_8beds_sen(flag, data_path, s0_all_matlab)

s0 = zeros(150,150,216); %zeros(301,301,431);
z=48; %95;

%
%data_path='/home/raid4/jnzhou/lmrecon_working/app/toshiba_real_data/01-29-2014/';
%
if flag

	for n=1:8

		s=touch(sprintf('%s/sensitivity_150x%d_%d.img',data_path, z,n));
		s=reshape(s, size(s0,1), size(s0,2), z);

		% put decay correction factor into sensitivity
	%	ff = sprintf('%s/bed_%d.mul_fac', data_path, n);
	%	m=touch(ff, 'float32', 1);
	%	decay_factor(n)=1/m;	
	%	nnn = decay_factor(n)/decay_factor(1);
		
		s0(:,:,[1:z]+24*(n-1)) = s0(:,:,[1:z]+24*(n-1)) + flipdim(s,3); % * nnn;

	end

else

	for n=1:8

		s=s0_all_matlab(:,:,:,n);
		% put decay correction factor into sensitivity
	%	ff = sprintf('%s/bed_%d.mul_fac', data_path, n);
	%	m=touch(ff, 'float32', 1);
	%	decay_factor(n)=1/m;	
	%	nnn = decay_factor(n)/decay_factor(1);
		
		s0(:,:,[1:z]+24*(n-1)) = s0(:,:,[1:z]+24*(n-1)) + flipdim(s,3); % * nnn;

	end

end