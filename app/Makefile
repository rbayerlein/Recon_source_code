CC=g++
INC=../inc
SRC_PATH=../src
SRC=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp

lmrecon:
	$(CC) -o lmrecon lmrecon.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm

lmsen_toshiba:
	$(CC) -o lmsen_toshiba lmsen_toshiba.cpp $(SRC) -I$(INC) -O3 -lm -fopenmp -DUSE_OMP -DUSE_SIDDON

lmdata_toshiba:
	$(CC) -o lmdata_toshiba lmdata_toshiba.cpp $(SRC) -O3 -I$(INC) 

prl2lm:
	$(CC) -o prl2lm myprl_decode.cpp -O3

nlm_flt:
	$(CC) -o nlm_flt nlm_flt.cpp -O3 -fopenmp -I$(INC)


