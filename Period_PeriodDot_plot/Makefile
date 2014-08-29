
DOWNLOADED = psrcat_pkg.tar
ATNF_DATABASE = psrcat_tar
DATA_FILE = name_p0_p1_AGE_BSURF.txt
PSRCAT_FILE_PATH = /home/greg/Neutron_star_modelling/Period_PeriodDot_plot/psrcat_tar/psrcat.db
SAVE_FIG = Period_PeriodDot.pdf

all:  $(DATA_FILE) $(ATNF_DATABASE) $(SAVE_FIG)

.PHONY: clean

$(ATNF_DATABASE): 
	wget http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz
	gunzip psrcat_pkg.tar.gz
	tar -xvf psrcat_pkg.tar
	cd psrcat_tar && ls && . makeit

$(DATA_FILE):
	export PSRCAT_FILE=${PSRCAT_FILE_PATH}; \
		./psrcat_tar/psrcat -c "name p0 p1 AGE BSURF" -l "p1<1" > ${DATA_FILE}

$(SAVE_FIG):make_plot.py
	python make_plot.py

clean:
	rm -f ${DOWNLOADED}
