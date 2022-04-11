

# SpydrPick
# : https://github.com/santeripuranen/SpydrPick

docker pull quay.io/biocontainers/spydrpick:1.2.0--h78a066a_0

docker run \
    -v ${PWD}:/project \
    -w /project \
    --cpuset-cpus 0-1 \
    -ti quay.io/biocontainers/spydrpick:1.2.0--h78a066a_0

outdir="05-coevolution/spydrpick"
mkdir -p $outdir
cat 02-processedData/NN_thai/* > $outdir/all.fasta

rootdir=$PWD
cd $outdir
SpydrPick \
    "05-coevolution/spydrpick/all.fasta" \
    --linear-genome \
    -t 1 \
    -o "all"
cd $rootdir

# plot results
Rscript --vanilla Scripts/coevolution/SpydrPick_plot.R


# fastcov
# : https://www.nature.com/articles/srep30425

# for each protein, concatenate with E
seqdir="02-processedData/AA_E_concat"
Efile="envelope protein E.fasta"
ls 02-processedData/AA > .tmp
for i in $(seq 1 $(sed -n '$=' .tmp )); do
    p="$(sed -n ${i}p .tmp )"
    if [[ "${p}" == "${Efile}" ]]; then
        echo 'skipping E itself...'
    else
        Rscript --vanilla Scripts/sequence_prep/concatSeqs.R \
            "${seqdir}/${p}" \
            "02-processedData/AA/${Efile}" \
            "02-processedData/AA/${p}"
    fi
done

# analyze coevolution
ls $seqdir > .tmp
for i in $(seq 1 $(sed -n '$=' .tmp )); do
    p="$(sed -n ${i}p .tmp )"
    outdir="05-coevolution/fastcov/${p%%.fasta}"
    mkdir -p $outdir
    ./bin/fastcov_linux \
        "${seqdir}/${p}" \
        -o "$outdir"
done

# plot results
Rscript --vanilla Scripts/coevolution/fastcov_plot.R