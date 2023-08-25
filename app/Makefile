CC=g++
INC=../inc
SRC_PATH=../src

SRC=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp


INC_CIP=../inc_cip
SRC_CIP_PATH=../src_cip
SRC_CIP=$(SRC_CIP_PATH)/petsys_cfg.cpp $(SRC_CIP_PATH)/petsys_cscanner.cpp $(SRC_CIP_PATH)/petsys_lmprj.cpp $(SRC_CIP_PATH)/petsys_ipsf.cpp $(SRC_CIP_PATH)/petsys_spmtx.cpp $(SRC_CIP_PATH)/petsys_reg.cpp $(SRC_CIP_PATH)/petsys_log.cpp $(SRC_CIP_PATH)/petsys_prj.cpp $(SRC_CIP_PATH)/petsys_timer.cpp $(SRC_CIP_PATH)/petsys_recon.cpp


SRC1=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj_any_noaddmul.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp


# SRC2=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp $(SRC_PATH)/petsys_ikernel.cpp



lmrecon_nontof:
	$(CC) -o lmrecon_nontof lmrecon.cpp $(SRC) -I$(INC) -O3 -DDEBUG -fopenmp -DUSE_OMP -lm   # -DUSE_TOF 

lmrecon_tof:
	$(CC) -o lmrecon_tof lmrecon.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -DUSE_TOF-fopenmp -DUSE_OMP -lm

#	$(CC) -o lmrecon lmrecon.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -lm


lmrecon_kernel_tof:
	$(CC) -o lmrecon_kernel_tof lmrecon_kernel.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm



lmrecon_tof_test:
	$(CC) -o lmrecon_tof_test lmrecon.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm

lmrecon_kernel_tof_test:
	$(CC) -o lmrecon_kernel_tof_test lmrecon_kernel.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm




lmrecon_cip_tof:
	$(CC) -o lmrecon_cip_tof lmrecon_cip.cpp $(SRC_CIP) -I$(INC_CIP) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm


lm_fp:
	$(CC) -o lm_fp lm_fp.cpp $(SRC1) -I$(INC) -O3 -DDEBUG -fopenmp -DUSE_OMP -lm  # -DUSE_NOADDMUL

lm_fp_tof:
	$(CC) -o lm_fp_tof lm_fp.cpp $(SRC1) -I$(INC) -O3 -DDEBUG -DUSE_TOF -fopenmp -DUSE_OMP -lm  # -DUSE_NOADDMU

lm_fp_exp:
	$(CC) -o lm_fp_exp lm_fp_exp.cpp $(SRC1) -I$(INC) -O3 -DDEBUG -fopenmp -DUSE_OMP -lm  # -DUSE_NOADDMUL


lmsen_toshiba:
	$(CC) -o lmsen_toshiba lmsen_toshiba.cpp $(SRC) -I$(INC) -O3 -lm -fopenmp -DUSE_OMP -DUSE_SIDDON

lmdata_toshiba:
	$(CC) -o lmdata_toshiba lmdata_toshiba.cpp $(SRC) -O3 -I$(INC) 

prl2lm:
	$(CC) -o prl2lm myprl_decode.cpp -O3

nlm_flt:
	$(CC) -o nlm_flt nlm_flt.cpp -O3 -fopenmp -I$(INC)


