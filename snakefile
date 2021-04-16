import sys
import yaml
import os
from Bio import SeqIO
from Bio.Seq import Seq

# Read config
configfile: "config.yaml"
_dir = config['outdir']
_qdir=_dir+'/query_genome'
_tdir=_dir+'/target_genome'

_threads=config['threads']

_query_genome=config['query_genome']
_target_genome=config['target_genome']

_query_genome_chr=config['query_chromosomes']
_target_genome_chr=config['target_chromosomes']


#read list of chromosomes to use
if not (_query_genome_chr == None or _query_genome_chr==""):
    with open(_query_genome_chr) as f:
        _query_genome_chr=f.read().splitlines()
else:
    # read from fasta
    print('Reading Query genome')
    _query_genome_chr=[]
    for record in SeqIO.parse(_query_genome, "fasta"):
        seqid=record.id
        _query_genome_chr.append(seqid)


if not (_target_genome_chr == None or _target_genome_chr==""):
    with open(_target_genome_chr) as f:
        _target_genome_chr=f.read().splitlines()
else:
    # read from fasta
    print('Reading Target genome')
    _target_genome_chr=[]
    for record in SeqIO.parse(_target_genome, "fasta"):
        seqid=record.id
        _target_genome_chr.append(seqid)



rule all:
    input:
        expand("{wd}/{query_genome_chr}.syn",wd=_dir,query_genome_chr=_query_genome_chr,qdir=_qdir,tdir=_tdir),
        expand("{wd}/final.syn",wd=_dir)

rule extract_chromosomes:
    group: "group1"
    input:
        qgenome={_query_genome},
        tgenome={_target_genome}

    output:
        expand("{qdir}/{query_genome_chr}.fasta",qdir=_qdir,query_genome_chr=_query_genome_chr),
        expand("{tdir}/target.fasta",tdir=_tdir)
    
    run:
        #create out dir
        if not os.path.exists(_dir):
            os.makedirs(_dir)
        if not os.path.exists(_qdir):
            os.makedirs(_qdir)
        if not os.path.exists(_tdir):
            os.makedirs(_tdir)

        # extract query chrmosomes in individual files
        total=0
        for record in SeqIO.parse(_query_genome, "fasta"):
            seqid=record.id
            if seqid in _query_genome_chr:
                outfasta=os.path.join(_qdir,seqid+'.fasta')
                output_handle = open(outfasta, "w")
                SeqIO.write(record, output_handle, "fasta")
                total=total+1
        print('Written {} query seqs'.format(total))


        total=0
        to_write=[]
        for record in SeqIO.parse(_target_genome, "fasta"):
            seqid=record.id
            if seqid in _target_genome_chr:
                to_write.append(record)
                total=total+1
        
        t_outfasta=os.path.join(_tdir,'target.fasta')
        SeqIO.write(to_write, t_outfasta, "fasta")
        print('Written {} target seqs'.format(total))


rule perform_synteny:
    input:
        q="{wd}/query_genome/{query_genome_chr}.fasta",
        t="{wd}/target_genome/target.fasta"
    output:
        "{wd}/{query_genome_chr}.syn"

    shell:
        """
        thisq=$(basename -s .fasta {input.q})
        echo "nucmer --mum -t {_threads} -p {wildcards.wd}/$thisq {input.q} {input.t}"
        nucmer --mum -t {_threads} -p {wildcards.wd}/$thisq {input.q} {input.t}
        show-coords -T -r -c -l {wildcards.wd}/$thisq.delta > {wildcards.wd}/$thisq.coords
        awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($4>$3) $8="+"; else $8="-"; print  $12,$1,$2,$13,$3,$4,$7,$8}}' {wildcards.wd}/$thisq.coords | awk 'BEGIN{{OFS="\t"}}{{if($5<$6){{start=$5;stop=$6}}else{{start=$6;stop=$5}};print $1,$2,$3,$4,start,stop,$7,$8}}' | sed '1,4d' > {wildcards.wd}/$thisq.syn

        """



rule merge_results:
    input:
        expand("{wd}/{query_genome_chr}.syn",wd=_dir,query_genome_chr=_query_genome_chr)
    output:
        "{wd}/final.syn"
    shell:
        """
        cat {input} > {wildcards.wd}/final.syn
        """



















