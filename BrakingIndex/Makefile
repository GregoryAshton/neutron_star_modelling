
DOWNLOADED = psrcat_pkg.tar
ATNF_DATABASE = psrcat_tar
DATA_FILE = ATNF_data_file.txt
PSRCAT_FILE_PATH = /home/greg/Neutron_star_modelling/BrakingIndex/psrcat_tar/psrcat.db

all:  $(ATNF_DATABASE) $(DATA_FILE) $(SAVE_FIG) 

.PHONY: clean

$(ATNF_DATABASE):
	wget http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz
	gunzip psrcat_pkg.tar.gz
	tar -xvf psrcat_pkg.tar
	cd psrcat_tar && ls && . makeit

$(DATA_FILE):
	export PSRCAT_FILE=${PSRCAT_FILE_PATH}; \
		./psrcat_tar/psrcat -c "name F0 F1 F2 BINARY Type Age" -nonumber -o short > ${DATA_FILE}

clean:
	rm -f ${DOWNLOADED}
