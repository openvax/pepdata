# download T-cell response assay results
echo
echo "Downloading IEDB T-cell response assay data"
echo 
wget -nc http://www.iedb.org/doc/tcell_compact.zip && \
	unzip tcell_compact.zip && \
	rm tcell_compact.zip && \
	mv tcell_compact.csv data

# download MHC binding assay results 
echo
echo "Downloading IEDB MHC Binding assay data"
echo 
wget -nc http://www.iedb.org/doc/elution_compact.zip && \
	unzip elution_compact.zip \
	&& rm elution_compact.zip \
	&& mv elution_compact.csv data

echo 
echo "Download Dana Farber tumor antigen list"
echo 
python df_tumor_antigens.py 
