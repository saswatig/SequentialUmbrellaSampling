	SUBROUTINE ENERGY_SUB(ICYCLE,ENERGY_SW,EN,PRESSURE,V_HCHI,V_WCA)
	IMPLICIT NONE
	
	INCLUDE 'INITIAL_PARAMETERS.inc'
	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'ENERGY_HARMONIC.inc'
	INCLUDE 'ENERGY_WCA.inc'
	INCLUDE 'MATRIX_Y.inc'
	INCLUDE 'MATRIX_X.inc'
	INCLUDE 'ENERGY_HCHI.inc'

	INTEGER ICYCLE,ENERGY_SW
	DOUBLE PRECISION EN,PRESSURE,V_HCHI

	DOUBLE PRECISION V_HARMONIC,P_HARMONIC
	DOUBLE PRECISION V_WCA,P_WCA
	DOUBLE PRECISION P_HCHI


	CALL ENERGY_HARMONIC_SUB(ENERGY_SW,V_HARMONIC,P_HARMONIC)
	CALL ENERGY_WCA_SUB(ENERGY_SW,V_WCA,P_WCA)
	CALL ENERGY_HCHI_SUB(ENERGY_SW,V_HCHI,P_HCHI)

	EN=V_HARMONIC+V_WCA+h*V_HCHI
	PRESSURE=P_HARMONIC+P_WCA-h*P_HCHI+n/(beta*boxx*boxy)




	RETURN
	END

	
