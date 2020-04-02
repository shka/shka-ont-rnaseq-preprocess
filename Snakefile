import gzip
import re
import subprocess

from Bio import SeqIO

configfile: 'config.yml'

hub_root = 'hub/'+ config['reference']['genome']['id']
threads = config['general']['threads']

##

rule PrepareSpikeInSequences:
    output: 'ref/spikein.fa.gz'
    params: config['reference']['spikein']['commands']
    shell: '{params} > {output}'

##

config['reference']['genome']['file'] = \
    'ref/'+config['reference']['genome']['url'].split('/')[-1]

rule DownloadGenomeSequences:
    output: config['reference']['genome']['file']
    params: config['reference']['genome']['url']
    shell: 'wget -P ref {params}'

##

rule PrepareGenomeSequences:
    input: rules.DownloadGenomeSequences.output
    output: 'ref/genome.fa.gz'
    run:
        proc = subprocess.Popen('pigz -c > '+output[0],
                                shell=True, stdin=subprocess.PIPE, text=True)
        with gzip.open(input[0], 'rt') as infh:
            for record in SeqIO.parse(infh, 'fasta'):
                if re.search('^chr', record.id):
                    SeqIO.write(record, proc.stdin, 'fasta')
        proc.stdin.close()
        proc.wait()

##

config['reference']['annotation']['file'] = \
    'ref/'+config['reference']['annotation']['url'].split('/')[-1]

rule DownloadAnnotation:
    output: config['reference']['annotation']['file']
    params: config['reference']['annotation']['url']
    shell: 'wget -P ref {params}'

##

rule PrepareAnnotation:
    input: rules.DownloadAnnotation.output
    output: 'ref/annotation.gff.gz'
    shell: "unpigz -c {input} | grep -v '^#' | grep ^chr | pigz -c > {output}"
        
##

# TODO: Although 'samtools faidx' asks pbgzip, compare pigz vs pbgzip.
rule PrepareReferenceSequences:
    input:
        rules.PrepareGenomeSequences.output,
        rules.PrepareSpikeInSequences.output
    output: 'ref/ref.fa.gz'
    threads: threads
    shell: 'unpigz -c {input} | bgzip -@ {threads} -c > {output}'

##

rule CheckReferenceSequenceLengths:
    input: rules.PrepareReferenceSequences.output
    output: 'ref/ref.fa.gz.fai'
    shell: 'samtools faidx {input}'

##

# NOTE: To build index for human genome, it requires more than 16GB RAM.
rule BuildMinimap2Index:
    input: rules.PrepareReferenceSequences.output
    output: 'ref/ref.mmi'
    params: config['minimap2']['index_options']
    threads: threads
    shell: 'unpigz -c {input} | minimap2 -t {threads} {params} -d {output} -'

##

tmp = []
for file in os.listdir(config['input']['directory']):
    if re.search('.+\.(fastq|fq)\.gz', file):
        tmp.append(config['input']['directory']+'/'+file)
config['input']['files'] = tmp

rule PrepareRawSequences:
    input: config['input']['files']
    output: 'out/raw.fastq.gz'
    run:
        proc_out = subprocess.Popen('pigz -c > '+output[0],
                                 shell=True, stdin=subprocess.PIPE)
        for file in input:
            proc_in = subprocess.Popen('unpigz -c '+file, shell=True,
                                       stdout=proc_out.stdin)
            proc_in.wait()
        proc_out.stdin.close()
        proc_out.wait()

##

# NOTE: At PolishCluster, racon 1.4.13 does not allow secondary
# alignments in the input BAM.
rule Minimap2RawSequences:
    input:
        index = rules.BuildMinimap2Index.output,
        fastq = rules.PrepareRawSequences.output
    output:
        bam = 'out/raw.bam'
    params: config['minimap2']['mapping_options']
    threads: threads
    shell: """
minimap2 -t {threads} -ax splice {params} {input.index} {input.fastq} \
| samtools view -b -u -F 2304 \
| samtools sort -@ {threads} - -o {output}
"""

##

rule FilterAlignments:
    input: rules.Minimap2RawSequences.output
    output: 'out/filtered.bam'
    params: config['minimap2']['filtering_options']
    threads: threads
    shell: 'samtools view -@ {threads} -b {params} {input} > {output}'

##

rule IndexBAM:
    input: '{prefix}.bam'
    output: '{prefix}.bam.bai'
    threads: threads
    shell: 'samtools index -@ {threads} {input}'

##

rule BAM2GFF:
    input: '{prefix}.bam'
    output: '{prefix}.gff.gz'
    threads: threads
    shell: 'spliced_bam2gff -t {threads} -s -M {input} | pigz -c > {output}'

##

rule ClusterGFF:
    input: 'out/filtered.gff.gz'
    output:
        gff = 'out/filtered.clustered.gff.gz',
        tsv = 'out/cluster_memberships.tsv.gz',
        fifo = temp('tmp/PinfishClusterGFF')
    params:
        c = config['pinfish']['minimum_cluster_size'],
        d = config['pinfish']['exon_boundary_tolerance'],
        e = config['pinfish']['terminal_exon_boundary_tolerance'],
        p = config['pinfish']['minimum_isoform_percentage']
    threads: threads
    shell: """
mkfifo {output.fifo} && pigz -c {output.fifo} > {output.tsv} &
unpigz -c {input} \
| cluster_gff -t {threads} \
    -c {params.c} -d {params.d} -e {params.e} -p {params.p} -a {output.fifo} \
| pigz -c > {output.gff}
wait
"""

##

rule CollapsePartials:
    input: 'out/{prefix}.gff.gz'
    output:
        gff = 'out/{prefix}.collapsed.gff.gz',
        tmp = temp('tmp/{prefix}.gff')
    params:
        d = config['pinfish']['exon_boundary_tolerance'],
        e = config['pinfish']['terminal_exon_boundary_tolerance'],
        f = config['pinfish']['first_exon_boundary_tolerance']
    threads: threads
    shell: """
unpigz -c {input} > {output.tmp} \
&& collapse_partials -t {threads} \
    -d {params.d} -e {params.e} -f {params.f} {output.tmp} \
| pigz -c > {output.gff}
"""

##

rule PolishClusters:
    input:
        tsv = rules.ClusterGFF.output.tsv,
        bam = rules.FilterAlignments.output
    output:
        fasta = 'out/filtered.clustered.polished.fasta.gz',
        fifo1 = temp('tmp/polish_clusters1'),
        fifo2 = temp('tmp/polish_clusters2')
    params:
        c = config['pinfish']['minimum_cluster_size']
    threads: threads
    shell: """
mkfifo {output.fifo1} && pigz -c {output.fifo1} > {output.fasta} &
mkfifo {output.fifo2} && unpigz -c {input.tsv} > {output.fifo2} &
polish_clusters -t {threads} \
    -c {params.c} -a {output.fifo2} -o {output.fifo1} {input.bam}
wait
"""

##

rule Minimap2PolishedClusters:
    input:
        index = rules.BuildMinimap2Index.output,
        query = rules.PolishClusters.output.fasta
    output: 'out/filtered.clustered.polished.bam'
    params: config['minimap2']['mapping_options']
    threads: threads
    shell: """
minimap2 -t {threads} -ax splice {params} {input.index} {input.query} \
| samtools view -b -u -F 2304 \
| samtools sort -@ {threads} - -o {output}
"""

##

rule GFFRead:
    input:
        genome = rules.PrepareReferenceSequences.output,
        gff = 'out/{prefix}.gff.gz'
    output:
        fasta = 'out/{prefix}.fasta.gz',
        tmp1 = temp('tmp/{prefix}1'),
        tmp2 = temp('tmp/{prefix}2')
    shell: """
unpigz -c {input.genome} > {output.tmp1} &
unpigz -c {input.gff} > {output.tmp2} &
wait
tmpfile=`dirname {output.fasta}`/`basename {output.fasta} .gz`
gffread -g {output.tmp1} -w $tmpfile {output.tmp2} && pigz $tmpfile
"""

##

rule GFFCompare:
    input:
        annotation = rules.PrepareAnnotation.output,
        gff = 'out/filtered.clustered.polished.collapsed.gff.gz'
    output:
        stats = 'out/gffcompare.stats',
        tmp1 = temp('tmp/gffcompare1'),
        tmp2 = temp('tmp/gffcompare2')
    shell: """\
unpigz -c {input.annotation} > {output.tmp1} &
unpigz -c {input.gff} > {output.tmp2} &
wait
tmpfile=`dirname {output.stats}`/`basename {output.stats} .stats`
gffcompare -r {output.tmp1} -R -M -C -K -o $tmpfile {output.tmp2}
"""

##

rule report:
    input:
        bam = rules.Minimap2RawSequences.output,
        bin = 'bin/report.el',
        template = 'out/report.org'
    output: 'out/report.html'
    shell: 'emacs -l {input.bin} --batch --kill'

##

rule BAM4Hub:
    input:
        bam = 'out/{prefix}.bam',
        fai = rules.CheckReferenceSequenceLengths.output
    output: hub_root + '/{prefix}.bam'
    threads: threads
    shell: """
samtools view -@ {threads} -b {input.bam} \
    `cut -f 1 {input.fai} | grep chr | paste -s -d ' '` \ > {output}
"""

##

rule BigGenePredAS:
    output: 'ref/bigGenePred.as'
    shell: """
wget -P ref https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as
"""

##

rule GFF3ToBigGenePred:
    input:
        gff = 'out/{prefix}.gff.gz',
        bgpas = rules.BigGenePredAS.output,
        chrsize = rules.CheckReferenceSequenceLengths.output
    output:
        gff = temp('tmp/{prefix}.gff'),
        genePred = temp('tmp/{prefix}.genePred'),
        preBigGenePred = temp('tmp/{prefix}.bb.pre'),
        bigGenePred = '%s/{prefix}.bb' % hub_root
    shell: """
unpigz -c {input.gff} | grep -E '^(#|chr)' > {output.gff} \
&& gtfToGenePred {output.gff} {output.genePred} \
&& genePredToBigGenePred {output.genePred} {output.preBigGenePred} \
&& bedToBigBed -type=bed12+8 -tab -as={input.bgpas} \
    {output.preBigGenePred} {input.chrsize} {output.bigGenePred}
"""

##

rule TrackDB:
    input:
        bai1 = hub_root + '/filtered.bam.bai',
        bai2 = hub_root + '/filtered.clustered.polished.bam.bai',
        bb = hub_root + '/filtered.clustered.polished.collapsed.bb'
    output: hub_root + '/trackDb.txt'
    shell: """
echo "
track filtered_bam
shortLabel Filt.Aligns
longLabel Alignments of the raw sequences after the filtereing.
bigDataUrl `basename {input.bai1} .bai`
type bam

track filtered_clustered_polished_bam
shortLabel Polis.Aligns
longLabel Alignments of the polished sequences after the filtering and the clustering.
bigDataUrl `basename {input.bai2} .bai`
type bam

track filtered_clustered_polished_collapsed_gff
shortLabel Polis.Trans
longLabel Transcripts identified by the clustering, the polishing and the collapsing.
bigDataUrl `basename {input.bb}`
type bigGenePred
" > {output}
"""

##

rule Genomes:
    input: rules.TrackDB.output
    output: 'hub/genomes.txt'
    params: config['reference']['genome']['id']
    run:
        with open(output[0], 'wt') as fh:
            fh.write("genome %s\ntrackDb %s/trackDb.txt" % (params, params))

##

config['ucsc']['genomesFile'] = 'genomes.txt'

rule Hub:
    input: rules.Genomes.output
    output: 'hub/hub.txt'
    run:
        with open(output[0], 'wt') as fh:
            for key in config['ucsc']:
                fh.write("%s %s\n" % ( key, config['ucsc'][key] ))

##

rule all:
    input:
        'out/filtered.clustered.collapsed.gff.gz',
        'out/filtered.clustered.polished.collapsed.fasta.gz',
        rules.GFFCompare.output.stats,
        'out/report.html',
        rules.Hub.output

