
DOWNLOADED = psrcat_pkg.tar
ATNF_DATABASE = psrcat_tar
DATA_FILE = name_p0_p1.txt

all:  $(DATA_FILE)

.PHONY: clean

$(ATNF_DATABASE): 
	wget http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz
	gunzip psrcat_pkg.tar.gz
	tar -xvf psrcat_pkg.tar
	cd psrcat_tar && ls && . makeit

$(DATA_FILE):
	export PSRCAT_FILE=/home/greg/Neutron_star_modelling/Period_PeriodDot_plot/psrcat_tar/psrcat.db
	./psrcat_tar/psrcat -c "name p0 p1" -l "p1 < 1.0" > name_p0_p1.txt

clean:
	rm -f ${DOWNLOADED}
