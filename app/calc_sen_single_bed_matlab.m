function senimg = cal_sen_single_bed_matlab(folder, image_size, voxel_size, bed_id)

%image_size = [150, 150, 48];
%voxel_size = 4.0775 * [1 1 1];
scanner_1 = buildPET('toshiba'); %scanner_1 = scanner_1.setProjector('Linterp');
scanner_1 = scanner_1.setNumberOfProjectionsPerAngle(353);

a = touch(sprintf('%s/Grp0_Frm0_Sub%d_atten_norm.raw', folder, bed_id));
a = reshape(a, 353, 320, 48*48);
a = permute(a, [1 3 2]);
senimg = toshiba_bp_sino_nonTOF_matlab(scanner_1, 1./a, image_size, voxel_size, 353);
