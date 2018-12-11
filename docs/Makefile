# Makefile
.PHONY:	SRS

all:	SRS MG MIS PS SVnVP UVnVP SVnVR UVnVR

SRS:
	cd SRS && make && cd ..

MG:
	cd Design/MG && make && cd ../..

MIS:
	cd Design/MIS && make && cd ../..

PS:
	cd ProblemStatement && make && cd ..

SVnVP:
	cd VnVPlan/SystVnVPlan && make && cd ../..

UVnVP:
	cd VnVPlan/UnitVnVPlan && make && cd ../..

SVnVR:
	cd VnVReport/SystVnVReport && make && cd ../..

UVnVR:
	cd VnVReport/UnitVnVReport && make && cd ../..


clean:	cleanSRS cleanMG cleanMIS cleanPS cleanSVnVP cleanUVnVP cleanSVnVR cleanUVnVR

cleanPS:
	cd ProblemStatement && make clean && cd ..

cleanSRS:
	cd SRS && make clean && cd ..

cleanMG:
	cd Design/MG && make clean && cd ..

cleanMIS:
	cd Design/MIS && make clean && cd ..

cleanSVnVP:
	cd VnVPlan/SystVnVPlan && make clean && cd ..

cleanUVnVP:
	cd VnVPlan/UnitVnVPlan && make clean && cd ..

cleanSVnVR:
	cd VnVReport/SystVnVReport && make clean && cd ..

cleanUVnVR:
	cd VnVReport/UnitVnVReport && make clean && cd ..