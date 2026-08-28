*[BED]: a plain text interval format: chromosome, start, end, and optional name, score and strand. 0-based and half-open
*[bedGraph]: a plain text format carrying one value per interval: chromosome, start, end, value. 0-based and half-open
*[bedMethyl]: the per-site methylation table modkit writes: a BED interval per modified position, with the counts behind the fraction
*[GFF3]: the General Feature Format: one annotated feature per line with its attributes in the ninth column. 1-based and inclusive
*[VCF]: the Variant Call Format: one line per position where a sample differs from the reference. 1-based
*[SAM]: the text form of an alignment file, as samtools view writes it. 1-based
*[BAM]: the compressed binary form of SAM. karyon does not read it: pipe it through samtools view
*[CRAM]: a reference-compressed alignment format. karyon does not read it: pipe it through samtools view
*[PAF]: the Pairwise mApping Format minimap2 writes: one line per alignment between two sequences
*[FASTA]: sequences as plain text, each one preceded by a line beginning with a greater-than sign
*[cytoBand]: the UCSC table of cytogenetic bands: chromosome, start, end, band name and stain
*[CIGAR]: the string in an alignment that says which bases matched, were inserted, deleted or clipped
*[SVTYPE]: the VCF attribute naming the kind of structural variant: a deletion, a duplication, an inversion or a translocation
*[SJ.out.tab]: the splice junction table an aligner writes: one intron per line, with the reads that crossed it. 1-based and inclusive on the intron
*[Newick]: a phylogeny written as nested parentheses, with branch lengths after colons
*[NEXUS]: a block-structured file that can carry a phylogeny along with the data behind it
*[NHX]: New Hampshire eXtended: Newick with named attributes attached to each branch
*[MSA]: a multiple sequence alignment: several sequences padded with gaps so that homologous positions line up in columns
*[ORF]: an open reading frame: a stretch between a start codon and the next stop codon in the same frame
*[codon]: three consecutive bases, the unit a ribosome reads to choose one amino acid
*[pN/pS]: the ratio of amino-acid-changing to silent variation, read against one because that is where selection is neutral
*[dN/dS]: the ratio of amino-acid-changing to silent substitution rates, read against one
*[GC skew]: how far the two strands differ in their G and C content, which changes sign at a replication origin
*[Tajima]: Tajima's D, a summary of how the frequencies of variants depart from what neutral evolution would give
*[hemimethylation]: a site methylated on one strand and not on the other
*[ideogram]: a chromosome drawn as its cytogenetic bands, the striped shape from a karyotype
*[cladogram]: a tree drawn to show only who is related to whom, with branch lengths carrying no meaning
*[phylogram]: a tree whose branch lengths are drawn to scale, so distance along a branch means evolutionary change
*[tanglegram]: two trees drawn face to face with their shared tips joined, so a disagreement between them shows as a crossing
*[unrooted]: a tree drawn without committing to which end is the ancestor
*[dendrogram]: a tree drawn beside rows to show how they cluster
*[synteny]: stretches of two genomes that are each other's counterparts, drawn as ribbons between them
*[dotplot]: the same pairwise alignments drawn as points on two axes, one genome per axis
*[pileup]: the aligned reads themselves, stacked in rows over the reference
*[lollipop]: a point event drawn as a stem with a head, so the position stays exact while the mark stays visible
*[skyline]: an estimate of population size through time, drawn as the steps the estimate is made of
*[phylodynamics]: reading a population's history out of the shape of a phylogeny of dated samples
*[half-open]: an interval that includes its start and excludes its end, so 0 to 10 is ten bases and the next interval starts at 10
*[minimap2]: an aligner for long reads and whole genomes, which writes PAF
*[samtools]: the standard toolkit for alignment files, which is where SAM text comes from
*[modkit]: the tool that turns modified-base calls into a bedMethyl table
*[Gubbins]: a tool that finds recombinant stretches in a bacterial alignment and writes them as GFF3 with the taxa that carry them
*[CNVkit]: a copy-number caller whose segment tables karyon reads
*[ASCAT]: a copy-number caller that separates the two parental copies
*[InterProScan]: the tool that annotates protein domains, whose table karyon reads
*[WebAssembly]: a compilation target that runs in a browser, which is how the playground runs the real program on the page
*[MSRV]: minimum supported Rust version: the oldest compiler the crate is tested against
