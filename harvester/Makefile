# COSMIC Mutation Data; can download from https://grch37-cancer.sanger.ac.uk/cosmic/download
# Note: we want GrCh37 build for now since most databases use this build.
COSMIC_MUTATIONS=CosmicMutantExport.tsv.gz

# Create COSMIC lookup table suitable for annotating variants with chrom/start/end/ref/alt
cosmic_lookup_table.tsv:
	gzcat $(COSMIC_MUTATIONS) | tail -n +2 | python cosmic_lookup_table.py create_table > temp.tsv
	(head -n 1 temp.tsv ; tail -n +2 temp.tsv | sort | uniq) > cosmic_lookup_table.tsv
	rm -f temp.tsv

# Create cgi_mut_benchmarking.tsv file.
cgi_mut_benchmarking.tsv:
	grep MUT cgi_biomarkers_20170208.tsv | \
	egrep -v "oncogenic mutation|insertion|deletion|wildtype|promoters|_variant" | \
	cut -f4 | sed 's/\:/	/g' > cgi_mut_benchmarking.tsv

# Create oncokb_mut_benchmarking.tsv file.
oncokb_mut_benchmarking.tsv:
	cut -f1-2 allActionableVariants.txt | \
	egrep -v "del|ins|Fusion|Oncogenic|mutations|dup|Wildtype|Amplification|Truncating|splice" \
	> oncokb_mut_benchmarking.tsv

# Benchmark file.
benchmark.txt: cosmic_lookup_table.tsv oncokb_mut_benchmarking.tsv
	 python cosmic_lookup_table.py benchmark --lookup-file cosmic_lookup_table.tsv --benchmark-file oncokb_mut_benchmarking.tsv > benchmark.txt

clean:
	rm -f cosmic_lookup_table.tsv cgi_mut_benchmarking.tsv oncokb_mut_benchmarking.tsv benchmark.txt
