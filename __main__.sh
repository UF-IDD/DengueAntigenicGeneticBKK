
# Compute pairwise antigenic distances between Thai viruses
# from cartography distances
Rscript --vanilla Scripts/antigenicDist_prep/thai_map.R


# Quantify expected antigenic distance resulting from
# measurement uncertainties
Rscript --vanilla Scripts/antigenicDist_uncertainty/getTiterTableWithNoise.R
for i in $(seq 1 100 ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla Scripts/antigenicDist_uncertainty/inferCoords.R
done




        #   Fit substitution model for each individual protein
        #   ..................................................

# There are 10 non-overlapping protein products of dengue
# which are further processed into 14 proteins
for i in $(seq 1 14 ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla Scripts/fitEffects/nnls.R thai_map AA
done

# summarize the fits
Rscript --vanilla \
    Scripts/plotEffects/summarizeFit.R \
    03-fits/nnls/thai_map/AA/adj_none/aasub




        #   Envelope protein (E)
        #   .....................

# Extract envelope protein from AA data
Rscript --vanilla \
    Scripts/sequence_prep/aa_sampleProteins2length.R \
    "02-processedData/AA_Envelope" \
    --outLength 495 \
    --nSamples 1 \
    -ip 3 \
    --sortsite

# make predictions (100-fold)
for i in $(seq 1 100 ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla Scripts/plotEffects/predictDist.R \
        thai_map \
        03-fits/nnls/thai_map/AA/adj_none/aasub/envelope_protein_E \
        --fasta 02-processedData/AA_Envelope/1.fasta
done

# examine fits
Rscript -e "rmarkdown::render(
    'Scripts/plotEffects/plotEffects_E.R'
    , output_dir = 'Reports'
    , output_file = 'Reports/effects_E.html'
    , knit_root_dir = getwd()
)"





        #   Find candidate protein(s) that harbor signals
        #   beyond linkage with E
        #   .....................


#   Generate random sequences to the length of each protein
#   to see if the protein has more signal than the average of the polyprotein
#   .........................................................................

nSamples=20

# generate sequences
for proteinLength in $(grep -oE "[0-9]+$" 02-processedData/AA_lengths.txt ) ; do
    Rscript --vanilla \
        Scripts/sequence_prep/aa_random2length.R \
        "02-processedData/AA_random2length_${proteinLength}" \
        ${proteinLength} \
        ${nSamples}
done
   
# fit nnls with these sequences
for d in $(ls -d 02-processedData/AA_random2length_* ) ; do
for i in $(seq 1 $nSamples ); do
    export SLURM_ARRAY_TASK_ID=$i
    seqdir=$(basename $d )
    Rscript --vanilla Scripts/fitEffects/nnls.R thai_map ${seqdir}
done
done

# summarize RMSE of the fits
for fitdir in $(ls -d 03-fits/nnls/thai_map/AA_random2length_*/adj_none/aasub ) ; do
    Rscript --vanilla \
        Scripts/plotEffects/summarizeFit.R \
        ${fitdir}
done


#   Downsample E to the size of each protein
#   to assess whether the protein has excess/equivalent/less signal
#   compared to E
#   .............

nSamples=20

# generate subsampled E sequences
for proteinLength in $(grep -oE "[0-9]+$" 02-processedData/AA_lengths.txt ) ; do
    Rscript --vanilla \
        Scripts/sequence_prep/aa_sample2length.R \
        "02-processedData/AA_sample2length_${proteinLength}" \
        ${proteinLength} \
        ${nSamples} \
        "envelope protein E"
done

# fit nnls with these sequences
for seqdir in $(ls -d 02-processedData/AA_sample2length_* ); do
for iSample in $(seq 1 ${nSamples} ); do
    export SLURM_ARRAY_TASK_ID=1
    Rscript --vanilla Scripts/fitEffects/nnls.R thai_map $(basename ${seqdir} )/${iSample}
done
done

# summarize RMSE of the fits
for fitdir in $(ls -d 03-fits/nnls/thai_map/AA_sample2length_*/* | grep -E "AA_sample2length_[0-9]+" ) ; do
    Rscript --vanilla \
        Scripts/plotEffects/summarizeFit.R \
        ${fitdir}/adj_none/aasub
done




        #   Where are the sites that harbor these additional signals
        #   beyond linkage with E
        #   .....................


nSamples=300

for proteinLength in 30 60 ; do
    # Generate random 30 or 60AA downsamples of NS2A
    Rscript --vanilla \
        Scripts/sequence_prep/aa_sampleProteins2length.R \
        "02-processedData/AA_NS2A_${proteinLength}" \
        --outLength ${proteinLength} \
        --nSamples ${nSamples} \
        -ip 7
    # Estimate the effect sizes, adjusted for E
    for i in $(seq 1 $nSamples ); do
        export SLURM_ARRAY_TASK_ID=$i
        Rscript --vanilla \
            Scripts/fitEffects/nnls.R \
            thai_map \
            AA_NS2A_${proteinLength} \
            --adjfasta 02-processedData/AA_Envelope/1.fasta \
            --adjname adj_envelope
    done
    
    # Generate random 30 or 60AA
    # from sites other than E and NS2A
    Rscript --vanilla \
        Scripts/sequence_prep/aa_sampleProteins2length.R \
        "02-processedData/AA_polyNonE_${proteinLength}" \
        --outLength ${proteinLength} \
        --nSamples ${nSamples} \
        -ip 1 5 6 8 9 10 12 11 14
    # Estimate the effect sizes, adjusted for E
    for i in $(seq 1 $nSamples ); do
        export SLURM_ARRAY_TASK_ID=$i
        Rscript --vanilla \
            Scripts/fitEffects/nnls.R \
            thai_map \
            AA_polyNonE_${proteinLength} \
            --adjfasta 02-processedData/AA_Envelope/1.fasta \
            --adjname adj_envelope
    done
done


# Summarize the fits
for proteinLength in 30 60 ; do
    # NS2A sites
    Rscript --vanilla \
        Scripts/plotEffects/summarizeFit.R \
        03-fits/nnls/thai_map/AA_NS2A_${proteinLength}/adj_envelope/aasub

    # Other sites
    Rscript --vanilla \
        Scripts/plotEffects/summarizeFit.R \
        03-fits/nnls/thai_map/AA_polyNonE_${proteinLength}/adj_envelope/aasub
done


# Generate report: 30AA
Rscript -e "rmarkdown::render(
    'Scripts/plotEffects/plotAdjustedEffects.R'
    , output_dir = 'Reports'
    , output_file = 'Reports/NS2A_30_adjE.html'
    , knit_root_dir = getwd()
    , params = list(
        fitdir = '03-fits/nnls/thai_map/AA_NS2A_30/adj_envelope/aasub'
        , nulldir = '03-fits/nnls/thai_map/AA_polyNonE_30/adj_envelope/aasub'
        , aadir = '02-processedData/AA_NS2A_30/'
        , lengthAdj = 495
    )
)"

# Generate report: 60AA
Rscript -e "rmarkdown::render(
    'Scripts/plotEffects/plotAdjustedEffects.R'
    , output_dir = 'Reports'
    , output_file = 'Reports/NS2A_60_adjE.html'
    , knit_root_dir = getwd()
    , params = list(
        fitdir = '03-fits/nnls/thai_map/AA_NS2A_60/adj_envelope/aasub'
        , nulldir = '03-fits/nnls/thai_map/AA_polyNonE_60/adj_envelope/aasub'
        , aadir = '02-processedData/AA_NS2A_60/'
        , lengthAdj = 495
    )
)"

# Compare positions identified from sampling 30 vs 60 AA
# : revealed 62 sites embedding signals beyond linkage with E
Rscript -e "rmarkdown::render(
    'Scripts/plotEffects/compareAdjustedEffects.R'
    , output_dir = 'Reports'
    , output_file = 'Reports/NS2A_compare_adjE.html'
    , knit_root_dir = getwd()
    , params = list(
        effectdir = '04-plots/fits_adjusted'
    )
)"




        #   Are these 62 sites mapped to the antigenic signals by chance?
        #   .............................................................

positionFile="04-plots/compare_adjusted/positions_NS2A_beyondE.txt"

#   E + 62 NS2A sites
#   .................

# Generate subset of NS2A sites that likely harbors signals beyond E
Rscript --vanilla \
    Scripts/sequence_prep/aa_subsetPositions.R \
    "02-processedData/AA_NS2A_beyondE" \
    -ip 7 \
    --positionFile ${positionFile}
# Fit the effect sizes, and get the performance measures
export SLURM_ARRAY_TASK_ID=1
Rscript --vanilla \
    Scripts/fitEffects/nnls.R \
    thai_map \
    AA_NS2A_beyondE \
    --adjfasta 02-processedData/AA_Envelope/1.fasta \
    --adjname adj_envelope
# Summarize results
Rscript --vanilla \
    Scripts/plotEffects/summarizeFit.R \
    "03-fits/nnls/thai_map/AA_NS2A_beyondE/adj_envelope/aasub"


#   E + 62 NS2A sites (permuted)
#   .................

# Permute NS2A sites that likely harbors signals beyond E
# (within site permutation across viruses)
# to break the non-random associations with diversity kept the same.
Rscript --vanilla \
    Scripts/sequence_prep/aa_permute.R \
    "02-processedData/AA_NS2A_beyondE_permuted" \
    --aaver "AA_NS2A_beyondE" \
    -ip 1 \
    -n 100 \
    --seed 58900
for i in $(seq 1 100 ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla \
        Scripts/fitEffects/nnls.R \
        thai_map \
        AA_NS2A_beyondE_permuted \
        --adjfasta 02-processedData/AA_Envelope/1.fasta \
        --adjname adj_envelope
done
# Summarize results
Rscript --vanilla \
    Scripts/plotEffects/summarizeFit.R \
    "03-fits/nnls/thai_map/AA_NS2A_beyondE_permuted/adj_envelope/aasub"



#   E + 62 other sites      
#   ..................

# Generate random AA with the same number as the number of sites
# suspected to harbor signals beyond E
# from sites other than E and NS2A
proteinLength=$( wc -l ${positionFile} | awk '{print $1}' )
Rscript --vanilla \
    Scripts/sequence_prep/aa_sampleProteins2length.R \
    "02-processedData/AA_polyNonE_${proteinLength}" \
    --outLength ${proteinLength} \
    --nSamples ${nSamples} \
    -ip 1 5 6 8 9 10 12 11 14
for i in $(seq 1 $nSamples ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla \
        Scripts/fitEffects/nnls.R \
        thai_map \
        AA_polyNonE_${proteinLength} \
        --adjfasta 02-processedData/AA_Envelope/1.fasta \
        --adjname adj_envelope
done
# Summarize results
export NOT_SAVE_ESTIMATES="TRUE" # do not save d,p,v estimates
Rscript --vanilla \
    Scripts/plotEffects/summarizeFit.R \
    "03-fits/nnls/thai_map/AA_polyNonE_${proteinLength}/adj_envelope/aasub"



        #   Sumarize the effect sizes of these 62 NS2A sites
        #   (adjusted for E)
        #   .....................................................

# make predictions (100-fold): when using E + 62 NS2A
for i in $(seq 1 100 ); do
    export SLURM_ARRAY_TASK_ID=$i
    Rscript --vanilla Scripts/plotEffects/predictDist.R \
        thai_map \
        03-fits/nnls/thai_map/AA_NS2A_beyondE/adj_envelope/aasub/subset \
        --fasta 02-processedData/AA_Envelope/1.fasta \
        --fasta 02-processedData/AA_NS2A_beyondE/subset.fasta
done

# Generate report: 62AA
Rscript -e "rmarkdown::render(
    'Scripts/plotEffects/plotSignalExcess.R'
    , output_dir = 'Reports'
    , output_file = 'Reports/NS2A_beyondE.html'
    , knit_root_dir = getwd()
    , params = list(plotdir = '04-plots/beyondE')
)"
