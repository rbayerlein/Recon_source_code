image_size = [150, 150, 48];
voxel_size = 4.0775 * [1 1 1];
scanner_1 = buildPET('toshiba');
scanner_1 = scanner_1.setNumberOfProjectionsPerAngle(353);

clear senimg;
for n=1:8, n

	a = touch(sprintf('%s/Grp0_Frm0_Sub%d_atten_norm.raw', data_path, n));
	a = reshape(a, 353, 320, 48*48);
	a = permute(a, [1 3 2]);
	senimg(:,:,:,n) = toshiba_bp_sino_nonTOF_matlab(scanner_1, 1./a, image_size, voxel_size, 353);

end

s0 = zeros(150,150,216); 
z=48; 
for n=1:8

	s=senimg(:,:,:,n);
	s0(:,:,[1:z]+24*(n-1)) = s0(:,:,[1:z]+24*(n-1)) + flipdim(s,3);

end
