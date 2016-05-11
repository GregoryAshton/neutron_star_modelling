
DOWNLOADED = psrcat_pkg.tar
ATNF_DATABASE = psrcat_tar
DATA_FILE = ATNF_data_file.txt
PSRCAT_FILE_PATH = /home/greg/neutron_star_modelling/Illustrations/Period_PeriodDot_plot/psrcat_tar/psrcat.db
SAVE_FIG = Period_PeriodDot.pdf

all:  $(ATNF_DATABASE) $(DATA_FILE) $(SAVE_FIG) 

.PHONY: clean

$(ATNF_DATABASE):
	wget http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz
	gunzip psrcat_pkg.tar.gz
	tar -xvf psrcat_pkg.tar
	cd psrcat_tar && ls && . makeit

$(DATA_FILE):
	export PSRCAT_FILE=${PSRCAT_FILE_PATH}; \
		./psrcat_tar/psrcat -c "name p0 p1 F2 BINARY Type Age" -nonumber -o short > ${DATA_FILE}

$(SAVE_FIG):make_plot.py
	python make_plot.py ${DATA_FILE}

clean:
	rm -f ${DOWNLOADED}
