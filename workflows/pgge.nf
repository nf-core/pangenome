process graphAligner {
  publishDir "${params.outdir}/testPGGE", mode: "${params.publish_dir_mode}"
    
  container "ghcr.io/pangenome/pgge:202103241056060dcc4b"

  input:
    path(fasta)
    // TODO find out if this each will execute things in parallel, I FEAR NOT
    each path(graph)

  output:
    path("${fasta}${graph}.gaf")

  """
  GraphAligner -g $graph -f $fasta -a ${fasta}${graph}.gaf -x vg -t ${task.cpus}
  """
}

process peanut {
  publishDir "${params.outdir}/testPGGE", mode: "${params.publish_dir_mode}"
    
  container "ghcr.io/pangenome/pgge:202103241056060dcc4b"

  input:
    path(gaf)

  output:
    path("${gaf}.pgge")

  shell:
  '''
  cut -f 2,3,4,16 !{gaf} | sed s/id:f:// | LC_NUMERIC=de_DE.UTF-8 awk '{ len=$3-$2; tlen+=len; sum+=$4*len; } END { print sum / tlen }' | tr "\n" "\t" > !{gaf}.pgge && peanut -g !{gaf} >> !{gaf}.pgge
  '''
}

workflow PGGE {
    take:
    fasta // channel: [val(f), path(fasta)]
    graphs // flattened channel: [graph1, graph2, ....]

    main:
        // TODO the Nextflow given splitFasta function might be slow, maybe we have to evaluate for larger sequences
        // splitFasta automatically unzips files
        // TODO WE WILL HAVE TO WRITE OUR OWN SPLIT FASTA BECAUSE THE FILE NAMES SUCK HARD
        ch_split_fasta = fasta.collect{it[1]}.splitFasta(file: true)

        graphAligner(ch_split_fasta, graphs)
        peanut(graphAligner.out)

    // emit:
        // some_output = samtoolsFaidx.out // path: *.txt
}