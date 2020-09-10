imax *
jmax *
kmax *
---------------
shapes * ch_#NAME# workspace_#NAME#.root #NAME#:$PROCESS
---------------
bin bin_#NAME#
observation -1
------------------------------
bin			ch_#NAME#	ch_#NAME#
process		sig			bkg_#NAME#
process		0			1
rate		#SIGINT#			#DATAINT#
--------------------------------
signorm	lnN		1.1		-
p0_#NAME# param		#P0#		#P0E#
p1_#NAME# param		#P1#		#P1E#
p2_#NAME# param		#P2#		#P2E#
